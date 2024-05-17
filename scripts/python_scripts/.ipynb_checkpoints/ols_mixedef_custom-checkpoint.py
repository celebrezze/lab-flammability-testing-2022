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

import warnings
warnings.filterwarnings("ignore")

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

def scale_and_center(df, datcolsall, cols_no_change=['cycle', 'doy', 'docy', 'cluster', 'csday1']):
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df[datcolsall])
    scaled_df = pd.DataFrame(scaled_data, columns=datcolsall)
    
    for col in cols_no_change:
        scaled_df[col] = df[col]
    
    return(scaled_df)

def formula_all_2way_interactions(cols, report=1, y='fh'):
    '''
    Returns just the maximal model for 2-way interactions - one formula including all 2-way interaction terms. (no singletons are inclus=ded bc all are accounted for in the interaction terms)
    '''
    cols2 = list(combinations(cols, 2))
    form = y+' ~ '
    for colset in cols2:
        form = form +colset[0]+'*'+colset[1]+' + '
    form = form[:-3]
    if report==1:
        print(form)
    return(form)

def mixedeff_check(df, col, yvar):
    '''
    Takes 1 predictor `X` and 1 outcome `y`. Performs OLS y = B0*X + B1. Prints model summary.
    '''
    form = yvar+'~'+col
    # model = smf.mixedlm(form, data=df, groups=df["plant_id"])
    model = smf.mixedlm(form, data=df, groups=df["species"], re_formula='1', vc_formula={'C(species):C(plant_id)': '0 + C(plant_id)'})
    results = model.fit(reml=False)
    Xcoef = results.params[1]
    pval = results.pvalues[1]
    aic_mod = results.aic
    return([Xcoef, pval, aic_mod])

def compare_predictors_mixedeff(dfog, cols, yvar='fh'):
    '''tests each col independently then reports all together in 1 table'''
    pvals = []
    coefs = []
    aics = []
    # compile coefs and pvals for each X col
    for col in cols:
        res = mixedeff_check(dfog, col, yvar)
        pvals.append(res[0])
        coefs.append(res[1])
        aics.append(res[2])
    df = pd.DataFrame({
        'cols':cols,
        'aics':aics,
        'pvals':pvals,
        'coefs':coefs
    })
    df['top_mod']=df.aics>=max(df.aics)-2
    df = df.sort_values('aics', ascending=False).reset_index(drop=True)
    print(df)

def compare_predictors_interaction_singletons(df, cols, y='fh', thresh=2, printsumm=0):

    aics = []
    colpairs = []
    int_terms = []
    
    # list of cols pairs
    cols2 = list(combinations(cols, 2))
    
    for col2 in cols2:
        
        # compile formula
        int_term = col2[0]+'*'+col2[1]
        formi = y+' ~ '+int_term
            
        # run model
        if printsumm==1:
            print(formi)
        try:
            model = smf.mixedlm(formi, data=df, groups=df["species"], re_formula='1', vc_formula={'C(species):C(plant_id)': '0 + C(plant_id)'})
            results = model.fit(reml=False)
            aics.append(results.aic)
            colpairs.append(col2)
            int_terms.append(int_term)
        except Exception as e:
            print("ERROR: Formula model error:", formi)

    # get interaction terms to keep by AIC
    df = pd.DataFrame({
        'colpairs':colpairs,
        'intterms':int_terms,
        'aics':aics,
    })
    df['top_mod']=df.aics<=min(df.aics)+thresh
    df = df.sort_values('aics', ascending=False).reset_index(drop=True)
    df_top = df[df.top_mod==True]

    # return and report
    sigcols = set([j for i in df_top.colpairs for j in i])
    sig_interactions = [interaction for interaction in df_top.intterms]
    print('\nColumns present in sig. interaction terms:', sigcols)
    print('\nTotal Num. Cols : Num. Sig. Int. Cols; ', len(cols), ':', len(sigcols))
    return sig_interactions

# GENERATE LISTS OF FORMULAS TO COMPARE

# Method given list of cols to use
def all_formulas_2way_interactions_and_singletons(cols_start, y='fh', report=1):
    '''
    Given a list of possible variables: 
    Return a list of all possible formulas built from all possible combinations of 2 way interaction terms and singletons which are not present in the interaction terms, including no singletons.
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

# Method given list of interactions and list of all other singletons
def red_formulas_2way_interactions_and_singletons(cols_start, base_interactions, y='fh', report=1):
    '''
    Given a list of known interaction terms and all variables: 
    Return a list of all possible formulas built from A LIST OF all possible combinations of 2 way interaction terms and singletons which are not present in the interaction terms, including no singletons.
    '''
    formulas = []
    colslist = []
    
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


def one_interaction_any_singletons():
    pass



# AIC ITERATION
def AICscore_from_all_pos_2way_interactions(df, formulas, report=0, thresh=2, rand_eff="plant_id"):
    '''
    Takes a list of formulas and a dataframe and returns a results dataframe of AIC scores, formulas, and columns used.
    '''
    scores = []
    formulas_return = []
    for formula in formulas:
        if report==1:
            print(formula)
        try:
            # get model & fit
            model = smf.mixedlm(formula, data=df, groups=df["species"], re_formula='1', vc_formula={'C(species):C(plant_id)': '0 + C(plant_id)'})
            results = model.fit(reml=False)
            # store score and formula
            scores.append(results.aic)
            formulas_return.append(formula)
        except Exception as e:
            print("ERROR: Formula model error:", formula)
            # pass
    
    # report best scores and formula
    # print(len(scores), len(formulas))
    resdf = pd.DataFrame({
            'AICscore':scores,
            'Formula':formulas_return
        }).sort_values(by='AICscore').reset_index(drop=True)
    num_top_models = len(resdf[resdf.AICscore<=resdf.AICscore.min()+thresh])
    if report==1:
        for i,row in resdf[resdf.AICscore<=resdf.AICscore.min()+thresh].iterrows():
            print(round(row.AICscore,2), ': ', row.Formula)
    return (resdf, num_top_models)


# FULL ITERATOR PROGRAM
def AIC_iterator(flam, cols_use, Y_VAR='fh',
                 minnumsingle=1, maxnumsingle=2, minnumint=0, maxnumint=1,
                 thresh_int=2, printsummint=0, #compare_predictors_interaction_singletons
                 report_AIC=0, thresh_AIC=2, rand_eff_AIC="plant_id" #AICscore_from_all_pos_2way_interactions
                ):
    '''
    Accepts dataframe and list of columns - from that list of columns it makes all possible formulas with the specified range of interaction and single terms. For interactions it only includes interactions determined to be significant in single testing. It includes all possible single terms regardles of standalone significance.
    '''
    
    # INTERACTIONS
    # significant singleton interactions: y = b + m1x1 + m2x2 + m3x1x2
    sig_interactions = compare_predictors_interaction_singletons(flam, cols_use, y=Y_VAR, thresh=thresh_int, printsumm=printsummint)
    # get list of tuples for interaction terms
    int_tuple_list = [tuple(x.split('*')) for x in sig_interactions]
    print('\nSignificant Interactions:')
    for pair in int_tuple_list:
        print(pair)

    # LIST OF FORMULAS
    df = flam
    cols = cols_use
    dv = Y_VAR
    
    formulas = []
    cols_used = []

    # BASE INTERACTIONS
    # all possible combinations of interactions from min to max
    intscombos = [list(combinations(int_tuple_list, n)) for n in range(minnumint, maxnumint+1)]

    # build form base starting with interaction terms
    for intcomboset in intscombos:
        for intcombo in intcomboset:
            # track cols used in interactions to drop from singletons later
            colsusedint = []
            # no interactions
            if len(intcombo)==0:
                form = dv+' ~ '
            # one or more interactions
            else:
                form = dv+' ~ '
                i = 0
                for int_tup in intcombo:
                    x1,x2 = int_tup
                    colsusedint.append(x1)
                    colsusedint.append(x2)
                    if i == 0:
                        form += x1+'*'+x2
                        i+=1
                    else:
                        form += ' + '+x1+'*'+x2
                # append base formula with only interactions
                formulas.append(form)
            # store list of cols used in interaction terms
            colsusedintset = set(colsusedint)
    
            # BEGIN ADDING SINGLES
            # create a copy of singletons list
            cols_wkg = cols.copy()
    
            # drop interactions terms from singletons list
            for xi in colsusedintset:
                cols_wkg.remove(xi)
            #print(colsusedintset)
            #print(cols_wkg)
    
            # generate list of all possible combos of singletons, from 1 to as many as there are
            singles_combos = [list(combinations(cols_wkg, n)) for n in range(minnumsingle, maxnumsingle+1)]
        
            # iterate over combo set (ie 1 poss singleton, 2 poss singletons, ... etc)
            for comboset in singles_combos:
                # for each combo in the combo set
                for combo in comboset:
                    # generate formula
                    form_wsingles = form
                    for single in combo:
                        if form_wsingles==dv+' ~ ':
                            form_wsingles+=single
                        else:
                            form_wsingles+=' + '+single
                    formulas.append(form_wsingles)
    
    print('\nNumber of formulas:', len(formulas))
            
    # AIC ITERATION
    resdf_fh, num_top_models = AICscore_from_all_pos_2way_interactions(df, formulas, report=report_AIC, thresh=thresh_AIC, rand_eff=rand_eff_AIC)

    top_resdf = resdf_fh[resdf_fh.AICscore<= resdf_fh.AICscore.min()+thresh_AIC]
    
    # report
    print('\n')
    #for idx,row in resdf_fh[0:num_top_models].iterrows():
    for idx,row in top_resdf.iterrows():
        formula = row.Formula
        print(formula)
    print('\n')
    #for idx,row in resdf_fh[0:num_top_models].iterrows():
    for idx,row in top_resdf.iterrows():
        formula = row.Formula
        model = smf.mixedlm(formula, data=df, groups=df["species"], re_formula='1', vc_formula={'C(species):C(plant_id)': '0 + C(plant_id)'})
        results = model.fit(reml=False)
        print(results.summary())
        plot_ols_coefficients(results)
        plt.show();
        # if 'species' in cols:
        #     cols.remove('species')
        # plot_resid(df, cols, results)
    
    



