import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
from itertools import chain
from sklearn.decomposition import PCA
# OLS
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from itertools import combinations, permutations
from sklearn.preprocessing import StandardScaler

# -----------------------------------------------------------------
# PLOTS
def corrplot(gdf, cols_of_interest):
    # correlation plot
    corgdf = gdf[cols_of_interest]
    cormtrx = corgdf.corr()
    # Create a mask to hide upper triangle
    mask = np.triu(np.ones_like(cormtrx, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(cormtrx, annot=True, cmap='coolwarm', fmt=".2f", mask=mask)
    plt.title('Correlation Matrix')
    plt.show();

def PCAplot(gdf, cols_of_interest, zoom=0.75, titletag='', dropnavals = 1):
    # Drop NA's?
    if dropnavals==1:
        gdf.dropna(subset=cols_of_interest, inplace=True)
    # Standardize the data
    scaler = StandardScaler()
    df = gdf[cols_of_interest]
    df_standardized = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)
    # Apply PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df_standardized)
    # Create a DataFrame with the PCA results
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    # Plot the PCA results with arrows representing the original features
    plt.figure(figsize=(10, 10))
    l = zoom
    plt.xlim(-1*l,l)
    plt.ylim(-1*l,l)
    plt.scatter(pca_df['PC1'], pca_df['PC2'], color='y', alpha=0.7)
    # Draw arrows for the original features
    for i, feature in enumerate(df.columns):
        plt.arrow(0, 0, pca.components_[0, i], pca.components_[1, i],
                  color='r', alpha=0.5, width=0.001, head_width=0.01, head_length=0.01)
        plt.text(pca.components_[0, i] * 1.2, pca.components_[1, i] * 1.2, 
                 feature, color='k', fontsize=14)
    plt.title('PCA Plot with Feature Arrows'+titletag)
    plt.xlabel('Principal Component 1 (PC1)')
    plt.ylabel('Principal Component 2 (PC2)')
    plt.grid(True)
    plt.show();

def regplot_smpl(x, y, o=3, ax=None):
    sns.regplot(x=x, y=y, ax=ax,
                scatter_kws={"color": 'r', "alpha": 0.1},
                line_kws={"color": 'k', "lw": 2},
                marker='+',  
                ci=95,
                order=o);  
    
def plot_resid(firesd1, cols, results):
    if len(cols)>9:
        print('ERROR: LENGTH COLS >9')
    residuals=results.resid
    fig, axs = plt.subplots(3, 3, figsize=(7, 6))  # 3 rows, 1 column
    axs_opt = list(permutations(range(3), 2))
    for i in range(len(cols)):
        a1,a2=axs_opt[i]
        regplot_smpl(x=firesd1[cols[i]], y=residuals, ax=axs[a1,a2])
        plt.tight_layout()
    plt.show();

def plot_ols_coefficients(results):
    """
    Plot OLS coefficients with error bars.
    Parameters:
    - results: OLS results object from statsmodels.
    Returns:
    - None (displays the plot).
    """
    coef = results.params
    conf_int = results.conf_int()
    fig, ax = plt.subplots(figsize=(10, 6))
    err = (conf_int.iloc[:, 1] - conf_int.iloc[:, 0]) / 2
    ax.bar(coef.index[1:], coef.values[1:], capsize=5, yerr = err.values[1:], 
           color='skyblue', edgecolor='black')
    ax.set_ylabel('Coefficient')
    ax.set_title('OLS Coefficients with Confidence Intervals')
    ax.set_xticklabels(coef.index[1:], rotation=45, ha="right")




# -----------------------------------------------------------------
# OLS

# 2-way interactions chain
def formula_all_2way_interactions(cols, report=1, y='fh'):
    '''
    Returns just the maximal model for 2-way interactions - one formula including all 2-way interaction terms. (no singletons are inclus=ded bc all are accounted for in the interaction terms)
    '''
    cols2= list(combinations(cols, 2))
    form = y+' ~ '
    for colset in cols2:
        form = form +colset[0]+'*'+colset[1]+' + '
    form = form[:-3]
    if report==1:
        print(form)
    return(form)

def all_formulas_2way_interactions_and_singletons(cols_start, y='fh', report=1):
    '''
    Return a list of all possible formulas built from all possible combinations of 2 way interaction terma and singletons which are not present in the interaction terms, including no singletons.
    '''
    formulas = []
    colslist = []
    # all options for 2-way interactions
    base_interactions = list(combinations(cols_start, 2))
    # for x in base_interactions:
    #     print(x)
    
    # all possible combinations of interactions from 1 to n possible interactions
    interaction_combos = [list(combinations(base_interactions, n)) for n in range(1, len(base_interactions)+1)]
            
    # assemble basis for formulas
    for combo_interact in interaction_combos:
        if report==1:
            print('************** interaction combo i *****************')
        # print(combo_interact)
        for x in combo_interact:
            x = list(x)
            # need to add in every possible combo of singles not already accounted for in the interaction
            cols_in_x = list(set(list(chain(*[list(t) for t in x]))))
            cols_not_in_x = [item for item in cols_start if item not in cols_in_x]
            singles_combos = [list(combinations(cols_not_in_x, n)) for n in range(0, len(cols_not_in_x)+1)]
            # add in singles
            for combo_sing in singles_combos:
                for s in combo_sing:
                    columns_used = cols_in_x
                    form_base = x + list(s)
                    if report==1:
                        print(form_base)
                    # assemble formulas
                    form = y + ' ~ '
                    for i in range(len(form_base)):
                        term = form_base[i]
                        if len(term)==2:
                            form = form + term[0] +'*'+ term[1] + ' + '
                            columns_used.append(term[0])
                            columns_used.append(term[0])
                        else:
                            form = form + term + ' + '
                            columns_used.append(term)
                    formulas.append(form[:-2])
                    colslist.append(list(set(columns_used)))
                    # report
                    if report==1:
                        print(form[:-2])
                        print(list(set(columns_used)))
    return formulas, colslist

def AICscore_from_all_pos_2way_interactions(df, formulas, form_cols, report=1, thresh=2, rand_eff="plant_id"):
    '''
    Takes a list of formulas and a dataframe and returns a results dataframe of AIC scores, formulas, and columns used.
    '''
    scores = []
    for formula in formulas:
        # get model & fit
        model = smf.mixedlm(formula, data=df, groups=df[rand_eff])
        results = model.fit(reml=False)
        # store score and formula
        scores.append(results.aic)
    
    # report best scores and formula
    resdf = pd.DataFrame({
            'AICscore':scores,
            'Formula':formulas,
            'form_cols':form_cols
        }).sort_values(by='AICscore').reset_index(drop=True)
    num_top_models = len(resdf[resdf.AICscore<resdf.AICscore.min()+thresh])
    if report==1:
        for i,row in resdf[resdf.AICscore<resdf.AICscore.min()+thresh].iterrows():
            print(round(row.AICscore,2), ': ', row.Formula)
    return (resdf, num_top_models)


def scale_and_center(df, datcolsall, cols_no_change=['cycle', 'doy', 'docy', 'cluster', 'csday1']):
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df[datcolsall])
    scaled_df = pd.DataFrame(scaled_data, columns=datcolsall)
    
    for col in cols_no_change:
        scaled_df[col] = df[col]
    
    return(scaled_df)

# 2-way interactions chain

# def all_poss_2way_interactions(cols, report=1):
#     '''
#     USES: formula_all_2way_interactions()
#     Takes list of cols to be tested, returns list of all possible 2-way interactions that can be made from those columns, from 1 interaction and up.
#     '''
#     # all possible combinations for interactions 2->all cols
#     combinations_2_or_more = []
#     for r in range(2, len(cols) + 1):
#         combinations_2_or_more.extend(combinations(cols, r))
#     if report==1:
#         print(combinations_2_or_more)
#     return(combinations_2_or_more)
# def AICselect_from_all_pos_2way_interactions(df, cols_start, report=1):
#     '''
#     USES: all_poss_2way_interactions()
#     Returns df of AIC scores & formulas for all possible 2-way interaction models
#     '''
#     # get all possible 2-way interaction combos
#     combos = all_poss_2way_interactions(cols_start, report=0)
#     # lists to store results
#     formulas =[]
#     form_cols = []
#     scores = []
#     # test all possible combos of 2-way interactions
#     for cols in combos:
#         # get formula
#         formula = formula_all_2way_interactions(cols, report=0)
#         # get model & fit
#         model = smf.ols(formula, data=df)
#         results = model.fit()
#         # store score and formula
#         scores.append(results.aic)
#         formulas.append(formula)
#         form_cols.append(cols)
#     # report best scores and formula
#     resdf = pd.DataFrame({
#         'AICscore':scores,
#         'Formula':formulas,
#         'form_cols':form_cols
#     }).sort_values(by='AICscore').reset_index(drop=True)
#     num_top_models = len(resdf[resdf.AICscore<resdf.AICscore.min()+2])
#     if report==1:
#         for i,row in resdf[resdf.AICscore<resdf.AICscore.min()+2].iterrows():
#             print(round(row.AICscore,2), ': ', row.Formula)
#     return (resdf, num_top_models)

def compare_corr_predictors_ols(dfog, cols, yvar='gr_rate'):
    pvals = []
    coefs = []
    # compile coefs and pvals for each X col
    for col in cols:
        res = linreg_check(dfog[col], dfog[yvar])
        coef = res[0]
        pval = res[1]
        pvals.append(pval)
        coefs.append(coef)
    df = pd.DataFrame({
        'cols':cols,
        'pvals':pvals,
        'coefs':coefs
    })
    df = df.sort_values('pvals')
    print(df.head())

def linreg_check(X, y):
    '''
    Takes 1 predictor `X` and 1 outcome `y`. Performs OLS y = B0*X + B1. Prints model summary.
    '''
    # Add a constant term to the independent variable matrix for statsmodels
    X = sm.add_constant(X)
    # Create and fit the OLS model
    model_sm = sm.OLS(y, X).fit()
    Xcoef = model_sm.params[1]
    pval = model_sm.pvalues[1]
    return([Xcoef, pval])












