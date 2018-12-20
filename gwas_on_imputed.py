from os.path import basename, join
from os import mkdir
from glob import glob
import bgen_reader
import numpy as np
import pandas as pd
import dask
from progressbar import ProgressBar as progbar
from progressbar import AdaptiveETA, Counter, Bar
from dask.diagnostics import ProgressBar


dtype = np.float32
samplefile = 'imputed/chr01.sample'
clean_index = pd.Index(np.loadtxt('data/clean_index.txt', dtype=str))

prefixes = ['healthspan']


sample = pd.read_csv(samplefile, skiprows=2, header=None, index_col=0, sep=' ')
mask = [True if elem in clean_index else False for elem in sample.index]
dfindexslice = sample.index[mask]

rho_df = pd.DataFrame(index=dfindexslice, dtype=dtype)
joshi_df = pd.DataFrame(index=dfindexslice, dtype=dtype)
Nd_df = pd.Series(dtype=dtype)
N_df = pd.Series(dtype=dtype)

def numpy_precalc(prefix):
    
    df_slice = pd.read_pickle(
        'data/df_slice_%s.pkl'%prefix).loc[dfindexslice,:].astype(dtype)
    betas = pd.read_pickle('data/betas_%s.pkl'%prefix).astype(dtype)
    M0, gamma = betas[['M0','gamma']].values
    N_ = 1.*len(df_slice)
    Nd_ = 1.*df_slice['event'].sum()

    df_slice['betax'] = df_slice[betas.index.drop(['M0','gamma'])].mul(betas, axis=1).sum(1)
    df_slice['z'] = np.exp(df_slice['betax'])
    df_slice['gt1'] = np.exp(df_slice['t1']*gamma)
    df_slice['gt2'] = np.exp(df_slice['t2']*gamma)
    df_slice['rho'] = np.exp(df_slice['betax'])*(np.exp(df_slice['t2']*gamma) - np.exp(df_slice['t1']*gamma))
    df_slice['rho'] /= df_slice['rho'].sum()
    df_slice['joshi'] = df_slice['event']/(df_slice['rho']*Nd_) - 1.

    rho_df.loc[:,prefix] = df_slice['rho'].values
    joshi_df.loc[:,prefix] = df_slice['joshi'].values
    Nd_df[prefix] = Nd_
    N_df[prefix] = N_


for prefix in prefixes:
    try:
        mkdir(join('data',prefix))
    except:
        pass
    numpy_precalc(prefix)


def make_dosage(f):
    bgen = bgen_reader.read_bgen(f, size=200, sample_file=samplefile, verbose=False)
    genotype = bgen['genotype'][:,mask,:].astype(dtype)
    dosage = genotype[:,:,1]+2.*genotype[:,:,2]
    return dosage, bgen


def calculate(dosage, rho, joshi, Nd, N):
    srho = (dosage*rho).sum(1)
    jorho = (joshi*rho).sum(1)[0]
    sjorho = ((dosage*joshi)*rho).sum(1)
    dsrho = dosage - srho.reshape((-1,1))
    dsrho2 = dsrho**2
    ds2rho = (dsrho2*rho).sum(1)

    beta1 = (sjorho - srho*jorho)/ds2rho
    sigma2 = (1./(Nd*ds2rho))**0.5

    return dask.array.stack([beta1, sigma2, dosage.sum(1)/(2.*N)], axis=1)


flst = np.sort(glob('imputed/chr*.bgen'))

pbar = progbar(widgets=[Bar(), Counter(), AdaptiveETA()])

for f in pbar(flst):

    dosage, bgen = make_dosage(f)
    variants = bgen['variants'].copy()

    resdfs = []

    for prefix in prefixes:

        rho = rho_df[prefix].values.reshape((1,-1))
        joshi = joshi_df[prefix].values.reshape((1,-1))
        Nd = Nd_df[prefix]
        N = N_df[prefix]

        resdfs.append(calculate(dosage, rho, joshi, Nd, N))


    with ProgressBar():
        x = dask.compute(*resdfs)

    for i, prefix in enumerate(prefixes):

        result = variants.copy()
        result['beta'] = x[i][:,0]
        result['sigma'] = x[i][:,1]
        result['eaf'] = x[i][:,2]
        result['N'] = N_df[prefix]
        result['Nd'] = Nd_df[prefix]
        outfile = join('data',prefix,
                       basename(f).replace('.bgen','.pkl'))
        result.set_index('id').to_pickle(outfile)