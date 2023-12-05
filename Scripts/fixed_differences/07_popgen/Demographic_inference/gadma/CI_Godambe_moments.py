#This is a script to simulate an SFS to produce confidence intervals using Godambe method. See the moments manual for more information.

#In this particular we are calculating CIs on the output from gadma_output_noMig10_th0 run.

import moments
import numpy as np

np.random.seed(seed=1000)

def model_func(params, ns):
	s1, t1, nu11, nu12 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	sts = moments.LinearSystem_1D.steady_state_1D(np.sum(ns))
	fs = moments.Spectrum(sts)
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	nu1_func = lambda t: (s1 * _Nanc_size) * (nu11 / (s1 * _Nanc_size)) ** (t / t1)
	nu2_func = lambda t: ((1 - s1) * _Nanc_size) * (nu12 / ((1 - s1) * _Nanc_size)) ** (t / t1)
	fs.integrate(tf=t1, Npop=lambda t: [nu1_func(t), nu2_func(t)], dt_fac=0.01)
	return fs

data = moments.Spectrum.from_file('/fixed_differences/07_popgen/Demographic_inference/gadma/CAT-SWE.sfs')
ns = data.sample_sizes

p0 = [0.8590291787561959, 0.18429451664439933, 1.2713137347673085, 5.427225544028091]
model = model_func(p0, ns)
ll_model = moments.Inference.ll_multinom(model, data)


fs_list = []

for i in range(100):
	fs_fixed = data.fixed_size_sample(np.rint(data.S()))
	fs_list.append(fs_fixed)

#print(fs_list)

uncerts_GIM = moments.Godambe.GIM_uncert(model_func, fs_list, p0, data)

print(uncerts_GIM)

# lower and upper CIs, in genetic units
lower = p0 - 1.96 * uncerts_GIM[:4]
upper = p0 + 1.96 * uncerts_GIM[:4]


print("95% CIs:")
print(f" s1 : {lower[0]:.1f} - {upper[0]:.1f}")
print(f" t1 : {lower[1]:.1f} - {upper[1]:.1f}")
print(f" nu11 : {lower[2]:.1f} - {upper[2]:.1f}")
print(f" nu12 : {lower[3]:.1f} - {upper[3]:.1f}")
