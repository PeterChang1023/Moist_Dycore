import numpy as np
import matplotlib.pyplot as plt

class ARStatisticalTest:
    def __init__(self, PC1_all, max_lag=100, num_samples=500):
        """
        Initialize the AR statistical test.

        Parameters:
        - PC1_all: (6, time) array, original PC1 time series for 6 experiments.
        - max_lag: int, maximum lag for autocorrelation calculation.
        - num_samples: int, number of Monte Carlo samples.
        """
        self.PC1_all = PC1_all
        self.max_lag = max_lag
        self.num_exp = PC1_all.shape[0]  # Number of experiments
        self.num_samples = num_samples

        # Store results
        self.auto_record = np.zeros((self.num_exp, self.max_lag + 1))
        self.residual = np.zeros_like(self.auto_record)
        self.AR_ensemble = np.zeros((self.num_exp, self.max_lag + 1, self.num_samples, PC1_all.shape[1]))
        self.auto_record_exp = np.zeros((self.num_exp, self.max_lag + 1, self.num_samples))

    def compute_autocorrelation(self):
        """
        Compute the autocorrelation of the original PC1 time series.
        """
        self.auto_record[:, 0] = 1  # ACF(0) = 1

        for lag in range(1, self.max_lag + 1):
            for exp in range(self.num_exp):
                self.auto_record[exp, lag] = np.corrcoef(
                    self.PC1_all[exp, :-lag], self.PC1_all[exp, lag:]
                )[0, 1]

                # Compute residual standard deviation
                self.residual[exp, lag] = np.std(
                    self.auto_record[exp, lag] * self.PC1_all[exp, :-lag] - self.PC1_all[exp, lag:]
                )

    def generate_AR_simulations(self):
        """
        Generate simulated time series using AR processes.
        """
        for exp in range(self.num_exp):
            for t in range(1, self.AR_ensemble.shape[3]):  # Time steps
                for lag in range(1, self.max_lag + 1):  # Lag values
                    self.AR_ensemble[exp, lag, :, t] = (
                        self.auto_record[exp, lag] * self.AR_ensemble[exp, lag, :, t-1] +
                        np.random.normal(0, scale=self.residual[exp, lag], size=self.num_samples)
                    )

    def compute_simulated_autocorrelation(self):
        """
        Compute the autocorrelation of the simulated AR time series.
        """
        for exp in range(self.num_exp):
            for lag in range(1, self.max_lag + 1):
                for ens in range(self.num_samples):
                    self.auto_record_exp[exp, lag, ens] = np.corrcoef(
                        self.AR_ensemble[exp, lag, ens, :-1], self.AR_ensemble[exp, lag, ens, 1:]
                    )[0, 1]

    def plot_results(self, input, exp_index=0):
        """
        input: usually is self.auto_record_exp
        exp_index=0: plot the dry exp, you can choose 0~5 mean 6 different exps
        Plot the boxplot of autocorrelation for different lags.
        """
        plt.figure(figsize=(10, 5))
        lags = np.arange(1, self.max_lag + 1)
        plt.boxplot(input[exp_index, 1:self.max_lag+1].T, positions=lags, widths=0.6, sym='')

        plt.xticks(np.linspace(0, self.max_lag, 11), labels=np.linspace(0, self.max_lag, 11, dtype=int))
        plt.xlabel("Lag (days)")
        plt.ylabel("Autocorrelation")
        plt.title(f"Boxplot of Autocorrelation for Different Lags (Exp {exp_index})")
        plt.grid()
        plt.show()

    def run_test(self):
        """
        Execute the full AR statistical test.
        """
        print("Computing original autocorrelation...")
        self.compute_autocorrelation()
        print("Generating AR simulations...")
        self.generate_AR_simulations()
        print("Computing simulated autocorrelation...")
        self.compute_simulated_autocorrelation()
        print("AR statistical test completed.")
        return self.auto_record_exp, self.AR_ensemble  # Return final simulated autocorrelation results

# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
from sklearn.decomposition import PCA

class EOF:
    """
    Calculating empirical orthogonal funcitons (EOFs)
    
    Parameters
    ----------
    dataset: tuple
        A tuple with elements are variables that you want to find their EOFs
        Variables must be array like, and must be standardized
        If given more than one dataset, combined EOF will be calculated
    
    n_components: int
        Number of modes that you need

    field: str, 1D or 2D, default = 2D
        The dimension of input variable arrays
    
    **svd_args: 
        Arguments for svd calculation in sklearn.decomposition.PCA
    
    About EOFs
    ----------
    The EOFs are vectors that represent the spatial distribution with largest temporal variation.
    In short, finding EOFs is equivalent to solving an eigenvalue problem of the variance matrix. The first eigen mode
    is EOF1, the second is EOF2, and so on.
    A variance matrix is done by multiplying the input variable array and its transpose, with temporal mean is zero.

    Note that
    ---------
    Original algorithm is developed by Kai-Chih Tseng: https://kuiper2000.github.io/
    """
    def __init__(
        self,
        dataset     : tuple,
        n_components: int,
        field       : str  = "2D",
        **svd_kwargs
    ):
        self.dataset      = dataset
        self.data_arr     = None
        self.n_components = n_components
        self.field        = field
        self.pca          = None
        self.EOF          = None
        self.PC           = None
        self.explained    = None
        self._svd         = svd_kwargs
    
    def _check_dimension(self):
        """
        If the dimensions of input variables are not consistent with self.field, raise ValueError

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        for sub in self.dataset:
            if (self.field == "2D" and np.ndim(sub) == 3) or (self.field == "1D" and np.ndim(sub) == 2): pass
            else:
                raise ValueError("The dimensions of input variables need to be consistent with input 'field'")

    def _single_subdataset_reshape_2D(self, subdataset: np.ndarray) -> np.ndarray:
        """
        Reshape input array with dimension (time, space) into (time*space)

        Parameters
        ----------
        subdataset: array
            The array of variable with dimension (time, space)
        
        Returns
        -------
        _subdataset_new: array
            The array of variable reshaped to dimension (time*space)
        """
        _subdataset_new = np.reshape(subdataset, (subdataset.shape[0], subdataset.shape[1]*subdataset.shape[2]))
        return _subdataset_new

    def _dataset_reshape_2D(self) -> tuple:
        """
        if there are more than two variables:
            Transfer input tuple with variable arrays into np.ndarray,
            and reshape it from dimension (var, time, space1, space2) into (time, var*space1*space2)
            Assign self.data_arr as the reshaped array
        else:
            Reshape the variable array into (time, space1*space2)
            Assign self.data_arr as the reshaped array

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if len(self.dataset) > 1:
            arr           = np.array(self.dataset)
            self.data_arr = np.reshape(np.transpose(arr, (1, 0, 2, 3)), (arr.shape[1], arr.shape[0]*arr.shape[2]*arr.shape[3]))
        else:
            self.data_arr = self._single_subdataset_reshape_2D(self.dataset[0])
    
    def _dataset_reshape_1D(self):
        """
        Same as _dataset_reshape_2D, but for 1-dimensional input variables

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if len(self.dataset) > 1:
            arr           = np.array(self.dataset)
            self.data_arr = np.reshape(np.transpose(arr, (1, 0, 2)), (arr.shape[1], arr.shape[0]*arr.shape[2]))
        else:
            self.data_arr = self.dataset[0]

    def _fit(self):
        """
        Create a PCA class and fit it with input data

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        pca_ = PCA(n_components = self.n_components, **self._svd)
        pca_.fit(self.data_arr)
        self.pca = pca_

    def _calc_EOF(self):
        """
        Calculate different EOF modes

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.EOF = self.pca.components_
    
    def _calc_PC(self):
        """
        Calculate PCs with input data and EOF modes

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        PC = np.dot(self.EOF, self.data_arr.T)
        self.PC = PC
    
    def _calc_explained(self):
        """
        Calculate the explainable ratio of each given EOF modes

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.explained = self.pca.explained_variance_ratio_

    def get(self):
        """
        Call _fit() _calc_EOF() _calc_PC _calc_explained() and calculate all of them

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self._check_dimension()
        if self.field == "1D":
            self._dataset_reshape_1D()
        else:
            self._dataset_reshape_2D()
        self._fit()
        self._calc_EOF()
        self._calc_PC()
        self._calc_explained()
##################################################################################################
import os
import gc
import h5py

class AtmosphericDiagnostics:
    """
    Perform general atmospheric calculations, including anomalies, 
    Eddy Momentum Flux (EMF), Eddy Heat Flux (EHF), and Eliassen-Palm (EP) flux.
    """

    def __init__(self, u, v, t, p, ps):
        """
        Initialize atmospheric diagnostics.

        Parameters:
        - u: ndarray, zonal wind (time, z, y, x)
        - v: ndarray, meridional wind (time, z, y, x)
        - t: ndarray, temperature (time, z, y, x)
        - p: ndarray, pressure levels (time, z, y, x)
        - ps: ndarray, surface pressure (time, y, x)
        """
        self.u = u
        self.v = v
        self.t = t
        self.p = p
        self.ps = ps

        self.time_dim = u.shape[0]
        self.z_dim = u.shape[1]
        self.y_dim = u.shape[2]
        self.x_dim = u.shape[3]

        self.lat = np.linspace(-90, 90, self.y_dim)
        self.lon = np.linspace(0, 360, self.x_dim)

    @staticmethod
    def anomaly(data, axis=-1):
        """Compute the anomaly of a given dataset by subtracting the mean along a specified axis."""
        return data - np.mean(data, axis=axis, keepdims=True)

    @staticmethod
    def EMF(u, v, axis=-1):
        """Compute the Eddy Momentum Flux (EMF) from wind components."""
        u_prime = AtmosphericDiagnostics.anomaly(u, axis=axis)
        v_prime = AtmosphericDiagnostics.anomaly(v, axis=axis)
        return u_prime * v_prime

    @staticmethod
    def EHF(t, p, v, axis=-1):
        """Compute the Eddy Heat Flux (EHF)."""
        Rd = 287.0
        Cp = 1004.0
        theta = t * (100000 / p) ** (Rd / Cp)

        v_prime = AtmosphericDiagnostics.anomaly(v, axis=axis)
        theta_prime = AtmosphericDiagnostics.anomaly(theta, axis=axis)

        b = (theta_prime / np.mean(theta, axis=axis, keepdims=True)) * 9.81
        b_prime = AtmosphericDiagnostics.anomaly(b, axis=axis)

        return v_prime * b_prime

    '''
    To EP flux...
    '''

    # @staticmethod
    # def compute_dtheta_dz(theta):
    #     """
    #     Compute the vertical gradient of potential temperature (dθ/dz).
        
    #     Parameters:
    #     - theta: ndarray, potential temperature (time, z, y, x)

    #     Returns:
    #     - dtheta_dz: ndarray, vertical gradient of theta (time, z, y, x)
    #     """
    #     return np.gradient(theta, axis=1)

    # @staticmethod
    # def EP_flux(f0, g, theta, p, v, emf, v_prime_b_prime):
    #     """
    #     Compute the Eliassen-Palm (EP) flux.

    #     Parameters:
    #     - f0: Coriolis parameter (scalar or array)
    #     - g: Gravitational acceleration (9.81 m/s²)
    #     - theta: ndarray, zonal mean potential temperature (time, z, y)
    #     - p: ndarray, pressure levels (time, z, y)
    #     - v: ndarray, meridional wind component (time, z, y, x)
    #     - emf: ndarray, Eddy Momentum Flux (time, z, y)
    #     - v_prime_b_prime: ndarray, meridional buoyancy flux (time, z, y)

    #     Returns:
    #     - F_j: Zonal component of the EP flux
    #     - F_k: Meridional component of the EP flux
    #     """
    #     dtheta_dz = np.gradient(theta, axis=1)  
    #     N_square = g / theta * dtheta_dz  

    #     F_j = emf.mean(axis=-1)  
    #     F_k = (f0[np.newaxis, np.newaxis, :] * v_prime_b_prime / N_square).mean(axis=-1)

    #     return F_j, F_k