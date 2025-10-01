import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, kendalltau, zscore
from sklearn.covariance import EllipticEnvelope
from sklearn.linear_model import RANSACRegressor
from scipy.optimize import curve_fit

# Generate sample data
#np.random.seed(42)
#x = np.random.rand(100)
#y = 2 * x + np.random.normal(0, 0.1, 100)
#
## Introduce outliers
#x[95:] += np.random.normal(0, 0.5, 5)
#y[95:] += np.random.normal(0, 0.5, 5)
#

def do(x,y,
       title_plot_pdf="Distribution",
       x_label_plot_pdf = 'Data',
       y_label_plot_pdf = 'No. of Occurences',
       x_label_distribution = 'Distribution x',
       y_label_distribution = 'Distribution y',
       title_plot_corr = "Correlation Methods and Outlier Detection",
       x_label_plot_corr = "X",
       y_label_plot_corr = "Y",
       figsize = [16,6],
       fig = None,
       namefig = None,
       methods2show=None, show=True):
    # Outlier detection using Z-score
    dim_data = x.shape[0]
    z_scores = np.abs(zscore(np.column_stack((x, y))))
    outliers_z = (z_scores >= 3).any(axis=1)
    x_z_filtered, y_z_filtered = x[~outliers_z], y[~outliers_z]
    print(f'Filtered {x_z_filtered.shape[0]}, original {dim_data}')
    perc_outliers_z = round(100 - 100*x_z_filtered.shape[0]/dim_data, 1)

    # Outlier detection using IQR
    q1_x, q3_x = np.percentile(x, [25, 75])
    iqr_x = q3_x - q1_x
    lower_x, upper_x = q1_x - 1.5 * iqr_x, q3_x + 1.5 * iqr_x
    q1_y, q3_y = np.percentile(y, [25, 75])
    iqr_y = q3_y - q1_y
    lower_y, upper_y = q1_y - 1.5 * iqr_y, q3_y + 1.5 * iqr_y
    outliers_iqr = (x < lower_x) | (x > upper_x) | (y < lower_y) | (y > upper_y)
    x_iqr_filtered, y_iqr_filtered = x[~outliers_iqr], y[~outliers_iqr]
    perc_outliers_iqr = round(100 - 100*x_iqr_filtered.shape[0]/dim_data, 1)

    # Outlier detection using EllipticEnvelope
    elliptic = EllipticEnvelope(contamination=0.1)
    X_elliptic = np.column_stack((x, y))
    #debug()
    try:
        elliptic.fit(X_elliptic)
        outliers_elliptic = elliptic.predict(X_elliptic) == -1
        x_elliptic_filtered, y_elliptic_filtered = x[~outliers_elliptic], y[~outliers_elliptic]
        perc_outliers_ellip = round(100 - 100*x_elliptic_filtered.shape[0]/dim_data, 1)
    except:
        outliers_elliptic=-1
        x_elliptic_filtered, y_elliptic_filtered = -1, -1
        perc_outliers_ellip =-1

    # RANSAC regression
    ransac = RANSACRegressor()
    try:
        ransac.fit(x.reshape(-1, 1), y)
        inliers_ransac = ransac.inlier_mask_
        x_ransac_filtered, y_ransac_filtered = x[inliers_ransac], y[inliers_ransac]
        perc_outliers_rsan = round(100 - 100*x_ransac_filtered.shape[0]/dim_data, 1)
    except:
        inliers_ransac=-1
        x_ransac_filtered, y_ransac_filtered = -1, -1
        perc_outliers_rsan = -1

    # Fit linear models
    def linear_model(x, a, b):
        return a * x + b

    fits = {}
    filtered_fits = {}
    filtered_datasets = {
        "Z-Score": (x_z_filtered, y_z_filtered),
        "IQR": (x_iqr_filtered, y_iqr_filtered),
        "EllipticEnvelope": (x_elliptic_filtered, y_elliptic_filtered),
        "RANSAC": (x_ransac_filtered, y_ransac_filtered)
    }
    methods = ["Pearson", "Spearman", "Kendall"]
    fit_statistics = {}
    filtered_fit_statistics = {}

    for method in methods:
        if method == "Pearson":
            params, cov = curve_fit(linear_model, x, y)
        else:
            params, cov = curve_fit(linear_model, np.sort(x), np.sort(y))
        fits[method] = params
        fit_statistics[method] = {"params": params, "std_err": np.sqrt(np.diag(cov))}
        
        filtered_fits[method] = {}
        filtered_fit_statistics[method] = {}
        
        for key, (x_f, y_f) in filtered_datasets.items():
            if isinstance(x_f,np.ndarray) and isinstance(y_f,np.ndarray):
                if method == "Pearson":
                    filtered_params, filtered_cov = curve_fit(linear_model, x_f, y_f)
                else:
                    filtered_params, filtered_cov = curve_fit(linear_model, np.sort(x_f), np.sort(y_f))
                filtered_fits[method][key] = filtered_params
                filtered_fit_statistics[method][key] = {"params": filtered_params, "std_err": np.sqrt(np.diag(filtered_cov))}
            else:
                filtered_params, filtered_cov = -1, -1 
                filtered_fits[method][key] = filtered_params
                filtered_fit_statistics[method][key] = {"params": filtered_params, "std_err": -1} 
           

    # Calculate correlations
    correlations = {}
    filtered_correlations = {}
    for method in methods:
        if method == "Pearson":
            corr, _ = pearsonr(x, y)
        elif method == "Spearman":
            corr, _ = spearmanr(x, y)
        elif method == "Kendall":
            corr, _ = kendalltau(x, y)
        correlations[method] = corr
        
        filtered_correlations[method] = {}
        for key, (x_f, y_f) in filtered_datasets.items():
            if isinstance(x_f,np.ndarray) and isinstance(y_f,np.ndarray):
                if method == "Pearson":
                    filtered_corr, _ = pearsonr(x_f, y_f)
                elif method == "Spearman":
                    filtered_corr, _ = spearmanr(x_f, y_f)
                elif method == "Kendall":
                    filtered_corr, _ = kendalltau(x_f, y_f)
                filtered_correlations[method][key] = filtered_corr
            else:
                filtered_correlations[method][key] = np.NaN


    # Plotting
    if show is True:
        if fig is None:
            fig, axes = plt.subplots(1, 2, figsize=figsize)
        else:
            axes_ori = fig.get_axes()
            # axes = [axes_ori[0], axes_ori[1]]
            axes = [axes_ori[2], axes_ori[3]]

        #try:
        if 1:
            # Left plot: Distribution with quartiles
            sns.histplot(x, ax=axes[0], kde=True, color='C0', alpha=0.5, label=x_label_distribution)
            sns.histplot(y, ax=axes[0], kde=True, color='C1', alpha=0.5, label=y_label_distribution)
            axes[0].axvline(lower_x, color='C0', ls='--')
            axes[0].axvline(upper_x, color='C0', ls='--')
            axes[0].axvline(lower_y, color='C1', ls='--')
            axes[0].axvline(upper_y, color='C1', ls='--')
            # For vlos
            #axes[0].axvline(0, color='w', ls='--')
            #amp4xlim = 3 
            #xlim4plots= [min([lower_x, lower_y])*amp4xlim, max([upper_x, upper_y])*amp4xlim]
            #axes[0].set_xlim(xlim4plots)
            # For vlos
            # For vturb
            #axes[0].axvline(0, color='w', ls='--')
            amp4xlim = 1.5 
            xlim4plots= [0, max([upper_x, upper_y])*amp4xlim]
            #xlim4plots= [0, 15]
            axes[0].set_xlim(xlim4plots)
            # For vturb
            axes[0].legend(fontsize=17)
            #axes[0].set_title("Distribution X and Y with Quartiles")
            axes[0].set_title(title_plot_pdf) 
            axes[0].set_xlabel(x_label_plot_pdf)
            axes[0].set_ylabel(y_label_plot_pdf)
            axes[0].margins(x=0, y=0)
            # Right plot: Scatter plot with original and filtered fits
            alpha = 0.5
            sns.scatterplot(x=x, y=y, label="Data Points", color="blue", alpha=0.3, ax=axes[1], s=30)
            #sns.kdeplot(x=x, y=y, fill=True, alpha=0.5, ax=axes[1], cmap='Reds') #, cmap='magma')
            #sns.scatterplot(x=x[outliers_z], y=y[outliers_z], label="Z-Score Outliers", color="purple", marker="X", ax=axes[1], alpha=alpha)
            #sns.scatterplot(x=x[outliers_iqr], y=y[outliers_iqr], label="IQR Outliers", color="brown", marker="D", ax=axes[1], alpha=alpha)
            #sns.scatterplot(x=x[outliers_elliptic], y=y[outliers_elliptic], label="EllipticEnvelope Outliers", color="cyan", marker="o", ax=axes[1], alpha=alpha)
            #sns.scatterplot(x=x[~inliers_ransac], y=y[~inliers_ransac], label="RANSACRegr. Outliers", color="magenta", marker="o", ax=axes[1], alpha=alpha, s=8)
            #    
            sns.scatterplot(x=x[outliers_z], y=y[outliers_z], label=f"Z-Score Outliers ({perc_outliers_z}%)", color="C1", marker=">", ax=axes[1], alpha=alpha, s=30)
            sns.scatterplot(x=x[outliers_iqr], y=y[outliers_iqr], label=f"IQR Outlier ({perc_outliers_iqr}%)", color="C3", marker="<", ax=axes[1], alpha=alpha, s=30)
            #sns.scatterplot(x=x[outliers_elliptic], y=y[outliers_elliptic], label=f"EllipticEnvelope Outliers ({perc_outliers_ellip}\%)", color="green", marker="^", ax=axes[1], alpha=alpha, s=30)
        #    sns.scatterplot(x=x[~inliers_ransac], y=y[~inliers_ransac], label=f"RANSACRegr. Outliers ({perc_outliers_rsan}%)", color="magenta", marker="o", ax=axes[1], alpha=alpha, s=8)


            colors = ["red", "green", "orange"]
            colors = ["C1", "C2", "magenta"]
            x_fit = np.linspace(min(x), max(x), 100)

            if methods2show is None:
                methods2show = methods.copy()
                print(f'WARNING: methods2show = {methods2show}')
            #for method, color in zip(methods, colors):
            count_color=0
            avoid_fiter = ['RANSAC', 'EllipticEnvelope']
            for method, qcolor in zip(methods2show, colors):
                y_fit = linear_model(x_fit, *fits[method])
                axes[1].plot(x_fit, y_fit, label=f"{method} Fit \n(r={correlations[method]:.2f}, a={fits[method][0]:.2f}, b={fits[method][1]:.2f})", color=f'C{count_color}', linestyle='-')
                count_color+=1
                count_color_2=0
                for key, (x_f, y_f) in filtered_datasets.items():
                    if isinstance(x_f,np.ndarray) and isinstance(y_f,np.ndarray):
                        print(key, (('RANSAC' in key) and ('EllipticEnvelope' in key)))
                        y_filtered_fit = linear_model(x_fit, *filtered_fits[method][key])
                        if (key in avoid_fiter) is False:
                            axes[1].plot(x_fit, y_filtered_fit, label=f"{method} {key} Fit \n(r={filtered_correlations[method][key]:.2f}, a={filtered_fits[method][key][0]:.2f}, b={filtered_fits[method][key][1]:.2f})", linestyle='--', color=colors[count_color_2])
                            count_color_2+=1
            axes[1].legend(fontsize=17, ncol=2)
            axes[1].set_title(title_plot_corr)
            axes[1].set_xlabel(x_label_plot_corr)
            axes[1].set_ylabel(y_label_plot_corr)
            xl = axes[1].get_xlim()
            axes[1].plot(xl,xl, ls='--', color='grey', alpha=.5)
            axes[1].margins(x=0, y=0)
            axes[1].set_xlim(xlim4plots)
            axes[1].set_ylim(xlim4plots)
            # For vturb
            #yl = axes[1].get_ylim()
            #axes[1].set_xlim([0, 15])
            #axes[1].set_ylim([0, 15])
            # For vturb
            plt.tight_layout()
        #except:
        #    print("ERROR ERROR ERROR: Exception error for plot!!!")

        if namefig is not None: plt.savefig(namefig)
        plt.show()


    # Return fit parameters, standard errors, and correlation factors
    return fits, filtered_fits, correlations, filtered_correlations, fit_statistics, filtered_fit_statistics, [perc_outliers_z, perc_outliers_iqr, perc_outliers_ellip, perc_outliers_rsan], [q1_x, q3_x, lower_x, upper_x, q1_y, q3_y, lower_y, upper_y]

