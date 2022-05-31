import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ez FIELD
Ez = S.Field.Field0.Ez(timesteps=300).getData()[0]
Validate("Ez field at iteration 300", Ez, 0.00002)

# TEST THAT Ubal_norm STAYS OK
uelm = S.Scalar.Uelm().getData(timestep=300)
Validate("Uelm", uelm, 1e-6 )

