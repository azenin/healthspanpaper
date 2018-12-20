from os.path import basename, join
from os import mkdir
import numpy as np
import pandas as pd
from pandas_plink import read_plink
from progressbar import ProgressBar as progbar
from progressbar import AdaptiveETA, Counter, Bar


def getbeta(inpsnpdata):

    notnansnp = (~np.isnan(inpsnpdata))
    snpdata = inpsnpdata[notnansnp]
    df = dfxs.iloc[notnansnp,:]

    delta, t1, t2, betax, z, gt1, gt2 = df.loc[:,['event','t1','t2','betax','z','gt1','gt2']].values.T

    g = z*(gt2-gt1)/gamma

    Nd = np.sum(delta)

    phi = M0*g - delta
    theta = np.copy(g)

    Ls = np.sum(snpdata*phi)
    Los = np.sum(snpdata*theta)
    Lss = M0*np.sum(snpdata**2*theta)
    Loo = (Nd/M0**2)

    commonpart = 1./(Lss-Los**2/Loo)
    beta = -Ls*commonpart
    sigma = np.sqrt(commonpart)

    return pd.Series(np.array([beta, sigma, Nd]),
                     index=['beta', 'sigma', 'Nd'])



prefix = 'healthspan'

try:
    mkdir(join('data',prefix))
except OSError:
    pass

df_slice = pd.read_pickle('data/df_slice_%s.pkl'%prefix).dropna(how='all')
betas = pd.read_pickle('data/betas_%s.pkl'%prefix)
M0, gamma = betas[['M0','gamma']].values
df_slice['betax'] = df_slice.loc[:,betas.index.drop(['M0','gamma'])].mul(betas, axis=1).sum(1)
df_slice['z'] = np.exp(df_slice['betax'])
df_slice['gt1'] = np.exp(df_slice['t1']*gamma)
df_slice['gt2'] = np.exp(df_slice['t2']*gamma)

for chrn in range(1,23):
    
#     chrn = '22'
    filename = 'calls/chr%02d' % chrn
    print filename

    (bim, fam, G) = read_plink(filename)

    mask = [True if elem in df_slice.index else False for elem in fam['fid']]
    xs = fam['fid'].values[mask]

    dfxs = df_slice.loc[xs,['event','t1','t2','betax','z','gt1','gt2']]

    g = G.rechunk((1000,-1))[:,mask].astype(float)

    pbar = progbar(widgets=[Bar(), Counter(), AdaptiveETA()])
    mass = []

    for block in pbar(list(g.blocks)):
        mass.append(pd.concat(map(getbeta, block.compute()), axis=1).T)


    result_table = pd.DataFrame(pd.concat(mass).values, index=bim['snp'], columns=['beta','sigma','Nd'])
    outfile = join('data', prefix, basename(filename)+'.pkl')
    result_table.to_pickle(outfile)