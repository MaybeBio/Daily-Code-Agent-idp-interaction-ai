import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import expon


def histogram_to_scatter(data, bins=50, density=True):
    """
    Convert histogram data to scatter plot coordinates.
    
    Parameters
    ----------
    data : array-like
        The data to histogram
    bins : int or array-like
        Number of bins or bin edges
    density : bool
        If True, normalize to probability density
        
    Returns
    -------
    x : array
        Bin centers
    y : array
        Bin heights (counts or densities)
    """
    counts, edges = np.histogram(data, bins=bins, density=density)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, counts


def exponential(x, rate, amplitude=1.0):
    """Exponential decay: f(x) = amplitude * rate * exp(-rate * x)"""
    return amplitude * rate * np.exp(-rate * x)


def power_law(x, alpha, amplitude=1.0):
    """Power law: f(x) = amplitude * x^(-alpha)"""
    return amplitude * x**(-alpha)


def fit_and_plot(data, fit_func=exponential, bins=50, 
                 figsize=(8, 5), log_scale=True, **kwargs):
    """
    Create scatter histogram and fit a function.
    
    Parameters
    ----------
    data : array-like
        Time series data (e.g., ON times)
    fit_func : callable
        Function to fit (default: exponential)
    bins : int
        Number of histogram bins
    figsize : tuple
        Figure size
    log_scale : bool
        Use log scale on y-axis
    **kwargs : dict
        Additional arguments passed to curve_fit (e.g., p0 for initial guess)
        
    Returns
    -------
    fig, ax : matplotlib figure and axis
    params : fitted parameters
    cov : covariance matrix
    """
    # Convert to scatter
    x, y = histogram_to_scatter(data, bins=bins, density=True)
    
    # Remove zero or negative values for fitting if using log scale
    mask = y > 0
    x_fit = x[mask]
    y_fit = y[mask]
    
    # Fit function
    try:
        params, cov = curve_fit(fit_func, x_fit, y_fit, **kwargs)
        fit_label = f"Fit: {fit_func.__name__}"
        
        # Generate smooth fit line
        x_smooth = np.linspace(x_fit.min(), x_fit.max(), 200)
        y_smooth = fit_func(x_smooth, *params)
        
    except Exception as e:
        print(f"Fitting failed: {e}")
        params, cov = None, None
        x_smooth, y_smooth = None, None
        fit_label = "Fit failed"
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.scatter(x, y, s=30, alpha=0.7, label='Data', color='C0', edgecolors='black', linewidth=0.5)
    
    if params is not None:
        ax.plot(x_smooth, y_smooth, 'r-', linewidth=2, label=fit_label)
    
    ax.set_xlabel('Time (frames)')
    ax.set_ylabel('Probability density')
    
    if log_scale:
        ax.set_yscale('log')
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Print fit parameters
    if params is not None:
        print(f"Fitted parameters: {params}")
        print(f"Covariance matrix:\n{cov}")
    
    return fig, ax, params, cov


# Example usage
if __name__ == "__main__":
    # Generate synthetic data (exponential distribution)
    np.random.seed(42)
    true_rate = 0.05
    synthetic_data = np.random.exponential(1/true_rate, size=1000)
    
    # Fit with exponential
    fig, ax, params, cov = fit_and_plot(
        synthetic_data, 
        fit_func=exponential,
        bins=30,
        p0=[0.1, 1.0]  # initial guess for [rate, amplitude]
    )
    
    ax.set_title('ON Time Distribution - Exponential Fit')
    plt.show()
    
    # Alternative: fit with power law
    fig2, ax2, params2, cov2 = fit_and_plot(
        synthetic_data, 
        fit_func=power_law,
        bins=30,
        p0=[2.0, 1.0]  # initial guess for [alpha, amplitude]
    )
    
    ax2.set_title('ON Time Distribution - Power Law Fit')
    plt.show()
