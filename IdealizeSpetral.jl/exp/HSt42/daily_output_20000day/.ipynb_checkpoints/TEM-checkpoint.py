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
import numpy as np
from scipy.interpolate import interp1d

class TEM:
    """
    Step 1. interpolation z -> z_new
    Step 2. Calculate v* and w*
    Step 3. Smooth v* and w* (there is singular point)
    Step 4. Calculate TEM streamfunction
    """

    ### Step 1. 
    def __init__(self, dataset_z):
        # # lat, lev
        # self.sigma_mean      = np.nanmean(p/ps, axis=(0,3))
        # self.sigma_onlyz     = np.nanmean(sigma_mean, axis=1)
        # self.y               = np.linspace(-90,90,64)
        # self.yy, self.sigma_mean2 = np.meshgrid(self.y,self.sigma_onlyz)

        # self.z = dataset_z
        print(dataset_z)

    def _interpolation(inut):
        z_new = np.zeros((((z.shape[1],22,64,128))))
        # theta_new = np.zeros((((z.shape[1],22,64,128))))
        
        for i in range(z.shape[1]):
            for j in range(64):
                  for k in range(128):
                    fe            = interp1d(np.linspace(0,20,20),z[0,i,:,j,k],  fill_value='extrapolate')
                    z_new[i,:,j,k] = fe(np.linspace(-1,21,22))
                    # fe            = interp1d(np.linspace(0,20,20),theta[i,:,j,k],  fill_value='extrapolate')
                    # theta_new[i,:,j,k] = fe(np.linspace(-1,21,22))
            if i % 100 == 0:
                print(f"int({i})")
        print("Step 1. done")

    def get(self):
        self._interpolation()
        
