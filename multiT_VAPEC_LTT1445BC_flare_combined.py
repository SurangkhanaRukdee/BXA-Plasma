## fit multi-T VAPEC used with LTT1445BC, call this in sherpa by %run -i *.py
## original file name: multiT_VAPEC_LTT1445BC_flare_combined_BXA4.py

import bxa.sherpa as bxa
from bxa.sherpa.background.pca import auto_background
from bxa.sherpa.background.models import ChandraBackground
from bxa.sherpa.background.fitters import SingleFitter

pc_to_cm = 3.086e18 
dist_cm = 6.9 * pc_to_cm

load_pha("ltt1445bc_flare_combined_src.pi")
set_stat('cstat')
#subtract()
#group_counts(20)
notice(0.5,2.5)

components = [xsvapec.p1, xsvapec.p2, xsvapec.p3, xsvapec.p4, xsvapec.p5, xsvapec.p6, xsvapec.p7, xsvapec.p8, xsvapec.p9, xsvapec.p10]

multitemp = xsvapec.p1 + xsvapec.p2 + xsvapec.p3 + xsvapec.p4 + xsvapec.p5 + xsvapec.p6 + xsvapec.p7 + xsvapec.p8 + xsvapec.p9 + xsvapec.p10

set_source(xsphabs.abs1 * multitemp)

id = 1
#fit background
# PCA background model: https://johannesbuchner.github.io/BXA/pca-background-models.html https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/example_pcabackground.py
convmodel = get_model(1)
bkg_model = auto_background(1)
set_full_model(convmodel + bkg_model*get_bkg_scale(1))

# empirical background model: https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/example_automatic_background_model.py
#fitter = SingleFitter(id, "27476_bc", ChandraBackground)
#try:
#	fitter.tryload()
#except IOError:
#	fitter.fit(plot=True)
#print('freezing background params')
#for p in get_bkg_model(id).pars: 
#	p.freeze()
#set_full_model(id, get_bkg_model(id) * get_bkg_scale(id) + get_response(id)(xsphabs.abs1 * xsvapec.p1))

abs1.nh = 1e-3
freeze(abs1)

for p in components[1:]:
	p.Ca = p1.Ca
	p.Ni = p1.Ni
	p.C = p1.C
	p.N = p1.N
	p.O = p1.O
	p.Ne = p1.Ne
	p.Mg = p1.Mg
	p.Al = p1.Al
	p.Si = p1.Si
	p.S = p1.S
	p.Ar = p1.Ar
	p.Fe = p1.Fe

# assign temperatures: 0.1, 0.3, 0.5, 0.7, ...
for p, kT in zip(components, np.logspace(-1, 1, 11)):
	print("component temperature: %.2f keV" % kT)
	p.kT = kT
	p.kT.min=0.99*kT
	p.kT.max=1.01*kT
	p.kT.freeze()


from sherpa.models.parameter import Parameter, UnaryOpParameter

# temperature of gaussian peak
kTpeak = Parameter('Tdist', 'kTpeak', 1, 0.1, 5, 0.1, 5)

# normalisation of gaussian peak
normpeak = Parameter('Tdist', 'normpeak', 0.001, 1e-6, 0.01, 1e-6, 0.01)

# width of gaussian peak
sigma = Parameter('Tdist', 'sigma', 0.2, 0, 2, 0, 2)

for i, p in enumerate(components):
	#p.norm = normpeak * exp(-(log(p.kT / kTpeak)/sigma)**2)
	p.norm = normpeak * UnaryOpParameter(-(UnaryOpParameter(p.kT / kTpeak, math.log10, 'log10')/sigma)**2, math.exp, 'exp')
	print("component norm: %.6f" % p.norm.val)

#UnaryOpParameter(p.kT, math.log, 'log')

print(get_model())
#-------------------------BXA-------------------------------

# where to store intermediate and final results? 
# this is the prefix used
outputfiles_basename = 'multiT_LTT1445BC_flare_combined_bxa4_bkgmodel'

# set model parameters ranges
p1.C.min = 0
p1.C.max = 1
p1.N.min = 0
p1.N.max = 1
p1.O.min = 0
p1.O.max = 1
p1.Ne.min = 0
p1.Ne.max = 1
p1.Mg.min = 0
p1.Mg.max = 1
p1.Al.min = 0 
p1.Al.max = 1
p1.Si.min = 0
p1.Si.max = 1
p1.S.min = 0
p1.S.max = 1
p1.Ar.min = 0
p1.Ar.max = 1
p1.Fe.min = 0
p1.Fe.max = 1


# list of parameters
parameters = [p1.C, p1.N, p1.O, p1.Ne, p1.Mg, p1.Al, p1.Si, \
p1.S, p1.Ar, p1.Fe, kTpeak, normpeak, sigma] 


# list of prior transforms	
priors = [
   bxa.create_uniform_prior_for(p1.C),
   bxa.create_uniform_prior_for(p1.N),
   bxa.create_uniform_prior_for(p1.O),
   bxa.create_uniform_prior_for(p1.Ne),
   bxa.create_uniform_prior_for(p1.Mg),
   bxa.create_uniform_prior_for(p1.Al),
   bxa.create_uniform_prior_for(p1.Si),
   bxa.create_uniform_prior_for(p1.S),
   bxa.create_uniform_prior_for(p1.Ar),
   bxa.create_uniform_prior_for(p1.Fe),
   bxa.create_loguniform_prior_for(kTpeak),
   bxa.create_loguniform_prior_for(normpeak),
   bxa.create_uniform_prior_for(sigma),
]
   
# make a single function:
priorfunction = bxa.create_prior_function(priors = priors)
print('running BXA ...')
solver = bxa.BXASolver(prior=priorfunction, parameters = parameters, 
	outputfiles_basename = outputfiles_basename)
results = solver.run(
	resume=True, # UltraNest supports resuming a crashed/aborted run
	frac_remain=0.5
	)

# plot inferred temperature distribution
import matplotlib.pyplot as plt
from ultranest.plot import PredictionBand
print ('Plotting kTdist')
plt.figure()
kT = np.linspace(components[0].kT.val, components[-1].kT.val, 400)
band = PredictionBand(kT)
for kTpeakval, normpeakval, sigmaval in solver.results['samples'][:,-3:]:
	band.add(normpeakval * np.exp(-(np.log10(kT / kTpeakval)/sigmaval)**2))
band.line(color='k')
band.shade(alpha=0.1, color='k')
plt.ylabel('Normalisation')
plt.xlabel('Temperature kT [keV]')
plt.xscale('log')
#plt.yscale('log')
plt.tight_layout()
plt.savefig(outputfiles_basename + '/plots/kTdist_1.pdf')
plt.close()

# plot inferred temperature distribution with emission measure
import matplotlib.pyplot as plt
from ultranest.plot import PredictionBand
print ('Plotting kTdist')
plt.figure()
kT = np.linspace(components[0].kT.val, components[-1].kT.val, 400)
band = PredictionBand(kT)
for kTpeakval, normpeakval, sigmaval in solver.results['samples'][:,-3:]:
	band.add(normpeakval / 1e-14 * (4 * np.pi * dist_cm**2) * np.exp(-(np.log10(kT / kTpeakval)/sigmaval)**2))
band.line(color='k')
band.shade(alpha=0.1, color='k')
plt.ylabel('Differential Emission Measure [cm$^{-5}$/keV]')
plt.xlabel('Temperature kT [keV]')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig(outputfiles_basename + '/plots/kTdist_2.pdf')
plt.close()

#exit()
print ('Calculating flux')
#--------------------------------Flux
import tqdm
prediction = []

for posterior_sample in tqdm.tqdm(solver.results['samples']):

	# load model (see above)
	for p, v in zip(solver.parameters, posterior_sample):
        	p.val = v


	# get plot data (see above)
	flux = calc_energy_flux(0.3,10.0, model=multitemp)
	#print(flux)


	# save plot data from this realisation:
	prediction.append(flux)

from uncertainties import ufloat
import uncertainties.unumpy as unumpy

flux = ufloat(np.median(prediction), np.std(prediction))
D = 6.9 * 3.1e+18
Lx = (4*numpy.pi*(D**2)*flux)
LogLx = unumpy.log10(4*numpy.pi*(D**2)*flux)

print ('Write down results')
#average flux from prediction list over the posteriors
np.savetxt(outputfiles_basename + '/flux_calc.txt', [[
	np.median(prediction),
	np.std(prediction),
	np.median(prediction)-np.quantile(prediction,0.158), #error bar of -1 sigma
	np.quantile(prediction,1-0.158)-np.median(prediction), #error bar of +1 sigma
	Lx,
	LogLx,
	]], delimiter=",",fmt='%e,%e,%e,%e,%s,%s', header='median, std, -1sigma, +1sigma, Lx, LogLx')

solver.set_best_fit()
print('cstat', calc_stat())
set_source(xsphabs.abs1 * multitemp)
set_stat('chi2gehrels')
group_counts(20)
ignore(None, 0.5)
ignore(2.5, None)

# need to have 3 additional free parameters to account for kTpeak, normalisation and width of the kT distribution		223	
p1.kT.thaw()		
p2.kT.thaw()		
p3.kT.thaw()		
# used the posterior results for the fit
p1.Fe = solver.results['posterior']['mean'][9]
p1.Ca = p1.Fe		
p1.Ni = p1.Fe		
p1.C = solver.results['posterior']['mean'][0]		
p1.N = solver.results['posterior']['mean'][1]	
p1.O = solver.results['posterior']['mean'][2]	
p1.Ne = solver.results['posterior']['mean'][3]		
p1.Mg = solver.results['posterior']['mean'][4]		
p1.Al = solver.results['posterior']['mean'][5]	
p1.Si = solver.results['posterior']['mean'][6]		
p1.S = solver.results['posterior']['mean'][7]		
p1.Ar = solver.results['posterior']['mean'][8]		

#subtract()
fit()
result = get_fit_results()
#result = get_stat_info()[0] # This get the best fit stat
set_analysis(1, 'energy', 'rate')
plt.figure(figsize=(5, 3))
plot_fit(color='black', clearwindow=True)
#result = get_stat_info()
#plt.xlim(0.2, 2.5)
plt.title('LTT1445BC flare combined')
Fe = solver.results['posterior']['mean'][9]

chi2 = '%.2f' % result.rstat
print(chi2)
#plt.text(2.0, 0.125,  '$\mathrm{kT_1}$ = %s keV' %kT1, size=10)
plt.text(0.98, 0.9, 'Red.$\chi^2$ = %s' % (chi2), size=10, transform=plt.gca().transAxes, va='top', ha='right')
plt.tight_layout()
plt.savefig('XXX.pdf', bbox_inches='tight')
plt.show()

plt.close()
plot_source()
plt.xlim(0.2, 2.5)
plt.savefig('XXX.pdf', bbox_inches='tight')
plt.show()
