import math
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys, os, shutil, argparse

### Description
# this script calculates persistence length of DNA from a DNAfold trajectory
  # file describing the thermal fluctuations of DNA over time.
# terminology: points (beads chosen for calculation after strand coarsening),
  # np (number of points), strand (chain of points, either beads or points,
  # along which the persistence length is calculated), nbd (number of beads).
# when giving simulation time, assumes dt=0.01.
# this script also estimates persistence length from E2E distribution between two
  # given points on the chain, as well as from fitting the log of tangent correlations.
# this version was created for use from the terminal.


def main():
	### for output in terminal
	print()

	### create argument parser
	parser = argparse.ArgumentParser()
	parser.add_argument('--file', 			default='output.dat',	help='trajectory file')
	parser.add_argument('--bases_skip',		default=0,	 			help='number of nucleotides to skip on each end')
	parser.add_argument('--steps_skip',		default=0,	 			help='number of recorded simulation steps to skip')
	parser.add_argument('--coarse_strand',	default=1,				help='coarse factor for nucleotides along strand')
	parser.add_argument('--coarse_time',	default=1,				help='coarse factor for time steps')

	### parse argument
	args = parser.parse_args()
	datFile = args.file
	nba_skip = int(args.bases_skip)
	steps_skip = int(args.steps_skip)
	coarse_strand = int(args.coarse_strand)
	coarse_time = int(args.coarse_time)

	### manually set input parameters
	legend_labels = ["Simulation"]								# list of legend labels
	MSD_bins = 30												# how many bins to divide the time into for MSD
	e2e_bins = 20												# how many bins to divide the time into for E2E

	### calculations
	if os.path.exists('figures'):
		shutil.rmtree('figures')
	os.mkdir('figures')
	lc_arr_nm,dot_avg_arr,ldot_sum_arr,e2e2_avg_arr,sim_time_ms,lp_est = calc_persistence(datFile,nba_skip,coarse_strand,steps_skip,coarse_time,MSD_bins,e2e_bins)

	### prepare figures
	lp_theory = 50
	fig1,axs = plt.subplots(2,2,figsize=(8,6))
	fig1.canvas.manager.set_window_title('Lp')
	fig1.tight_layout(pad=3)

	### calculate bond angle method
	lc_short = np.array([0])
	for j in range(1,len(lc_arr_nm)):
		if lc_arr_nm[j] < lp_theory/2:
			lc_short = np.hstack((lc_short,lc_arr_nm[j]))
	lcs_arr_nm = lc_short
	bam = np.array([0])
	for j in range(1,len(lc_short)):
		bam = np.hstack((bam,lc_arr_nm[j]/(1-dot_avg_arr[j])))
	bam_avg_arr = bam

	### plot simulation
	axs[0,0].semilogy(lc_arr_nm,dot_avg_arr,'--')
	axs[0,1].plot(lc_arr_nm,ldot_sum_arr,'--')
	axs[1,0].plot(lcs_arr_nm,bam_avg_arr,'--')
	axs[1,1].plot(lc_arr_nm,e2e2_avg_arr,'--')

	### prepare to plot log(correlation) theory
	legend_labels.append("Fit (P = " + str(round(lp_est)) + ")")
	legend_labels.append("P = " + str(lp_theory))

	### theory for dot
	fit_dot = np.exp(-lc_arr_nm/lp_est)
	theory_dot = np.exp(-lc_arr_nm/lp_theory)
	axs[0,0].semilogy(lc_arr_nm,fit_dot,'C0')
	axs[0,0].semilogy(lc_arr_nm,theory_dot,'k')
	axs[0,0].yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
	axs[0,0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

	### theory for sum of dots
	theory_sum = np.zeros((len(theory_dot),1))
	l0_avg = lc_arr_nm[1]-lc_arr_nm[0]
	for i in range(1,len(theory_dot)):
		theory_sum[i] = sum(theory_dot[1:(i+1)]) * l0_avg
	axs[0,1].plot(lc_arr_nm,theory_sum,color='k')
	axs[0,1].plot(lc_arr_nm,lc_arr_nm*0+lp_theory,color='k')

	### theory for l/(1-dot)
	axs[1,0].plot(lcs_arr_nm,np.zeros(len(lcs_arr_nm))+lp_theory,color='k')	

	### theory for end-to-end distance
	theory_e2e2 = 2*lc_arr_nm*lp_theory - 2*lp_theory**2*(1-np.exp(-lc_arr_nm/lp_theory))
	axs[1,1].plot(lc_arr_nm,theory_e2e2,'k')

	### plot formatting
	fig1.suptitle("Persistence Length from " + str(precision_round(sim_time_ms,5)) + " ms Simulation")
	axs[0,0].legend(legend_labels)
	axs[0,0].set(xlabel="Contour Length [nm]")
	axs[0,0].set(ylabel=r'$\langle cos\theta \rangle$')
	axs[0,1].set(xlabel="Contour Length [nm]")
	axs[0,1].set(ylabel=r'$\Sigma \ell \langle cos\theta \rangle$')
	axs[1,0].set(xlabel="Contour Length [nm]")
	axs[1,0].set(ylabel=r'$\ell / (1-\langle cos\theta \rangle$)')
	axs[1,1].set(xlabel="Contour Length [nm]")
	axs[1,1].set(ylabel=r'$\langle R^2 \rangle$')
	plt.savefig('figures/lp')

	### for output in terminal
	print()


def read_dnafold(fileName,nba_skip,coarse_strand,steps_skip,coarse_time):
	print("Loading DNAfold trajectory...")
	with open(fileName) as f:
		content = f.readlines()
		print("Parsing trajectory...")

		### bead information
		nbd_total = int(content[3].split()[0])
		nbd_scaf = 0
		while int(content[nbd_scaf+9].split()[1]) == 1:
			nbd_scaf += 1
		nbd_skip = int(np.floor(nba_skip/8))

		dbox = 2*int(content[5].split()[1])	
		np_scaf = int(np.floor((nbd_scaf-nbd_skip*2)/coarse_strand))
		steps_total = int(len(content)/(nbd_total+9))
		steps_use = int(np.floor((steps_total-steps_skip)/coarse_time))

		if np_scaf <=0:
			print("Too much cut off ends, using the entire chain")
			np_scaf = int(np.floor(nbd_scaf/coarse_strand))
			nbd_skip = 0

		dt_per_output = int(content[9+nbd_total+1].split()[0])
		sim_time_ms = (steps_total-steps_skip)*0.01*1E-6*dt_per_output*3300

		print("{:1.2e} steps in simulation".format(steps_total*dt_per_output))
		print("{:1.2e} steps recorded".format(steps_total))
		print("{:1.2e} steps in ensemble average".format(steps_use))

		v_tan = np.zeros((steps_use,np_scaf,3))
		points = np.zeros((steps_use,np_scaf,3))
		for i in range(steps_use):
			pbc = np.zeros(3)
			for j in range(np_scaf):
				line0 = content[(nbd_total+9)*(steps_skip+i*coarse_time)+9+j*coarse_strand+nbd_skip].split()
				line1 = content[(nbd_total+9)*(steps_skip+i*coarse_time)+9+j*coarse_strand+nbd_skip+1].split()
				if coarse_strand != 1:
					lineN = content[(nbd_total+9)*(steps_skip+i*coarse_time)+9+(j+1)*coarse_strand+nbd_skip].split()
				else:
					lineN = line1
				for k in range(3):
					points[i,j,k] = float(line0[k+2])*dbox - pbc[k]*dbox
					v_tan[i,j,k] = applyPBC( float(line1[k+2])*dbox - float(line0[k+2])*dbox ,dbox)
					pbc[k] = pbc[k] + int(round( (float(lineN[k+2])*dbox - float(line0[k+2])*dbox) /dbox))
			if coarse_strand == 1 and nbd_skip == 0:
				v_tan[i,j,:] = v_tan[i,j-1,:]
	return np_scaf,steps_use,v_tan,points,dbox,sim_time_ms


def calc_dot_avg(v_tan,points,np_scaf,steps_use):
	print("Calculating correlations...")
	dot_avg_arr = np.array([1])
	ldot_sum_arr = np.array([0])
	e2e2_avg_arr = np.array([0])
	l0_avg = calc_l0_avg(points)
	for di in range(1,np_scaf):
		dot_sum = 0
		e2e2_sum = 0
		for iL in range(np_scaf-di):
			for j in range(steps_use):
				iR = iL + di
				vL = [v_tan[j,iL,0],v_tan[j,iL,1],v_tan[j,iL,2]]
				vR = [v_tan[j,iR,0],v_tan[j,iR,1],v_tan[j,iR,2]]
				pL = np.array([points[j,iL,0],points[j,iL,1],points[j,iL,2]])
				pR = np.array([points[j,iR,0],points[j,iR,1],points[j,iR,2]])
				dot_sum += calc_dot(vL,vR)
				e2e2_sum += (np.linalg.norm(pR-pL))**2
		dot_avg_arr = np.hstack((dot_avg_arr,dot_sum/(steps_use*(np_scaf-di))))
		ldot_sum_arr = np.hstack((ldot_sum_arr,ldot_sum_arr[di-1]+dot_avg_arr[di]*l0_avg))
		e2e2_avg_arr = np.hstack((e2e2_avg_arr,e2e2_sum/(steps_use*(np_scaf-di))))
	return dot_avg_arr,ldot_sum_arr,e2e2_avg_arr,l0_avg


def lp_from_correlation(lc_arr_nm,dot_avg_arr):
	flag = False
	for i in range(len(dot_avg_arr)):
		if dot_avg_arr[i] <= 0:
			flag = True
			print("Failed to estimate persistence length from fitting log(correlation)")
			Lp_guess = 1			
			break
	if flag == False:
		slope = stats.linregress(lc_arr_nm,np.log(dot_avg_arr))[0]
		Lp_guess = -1/slope
	return Lp_guess


def e2e_dist(points,l0_avg,e2e_first,e2e_last,e2e_bins):
	steps_use = int(points.shape[0])
	np_scaf = int(points.shape[1])
	if e2e_last > np_scaf-1:
		print("Requested e2e_last too large, using last bead in chain.")
		e2e_last = np_scaf-1
	e2e_arr = np.zeros(steps_use)
	for i in range(steps_use):
		pL = np.array([points[i,e2e_first,0],points[i,e2e_first,1],points[i,e2e_first,2]])
		pR = np.array([points[i,e2e_last,0],points[i,e2e_last,1],points[i,e2e_last,2]])
		e2e_arr[i] = (np.linalg.norm(pR-pL))
	Lc = l0_avg*(e2e_last-e2e_first)
	e2e_avg = math.sqrt(np.mean(e2e_arr**2))

	### estimate persistence length
	Lp_guess = 100
	threshold = 0.001
	max_iter = 1000
	error = 1
	counter = 0
	while abs(error) > threshold:
		e2e = np.sqrt(2*Lp_guess*Lc - 2*Lp_guess**2* (1-math.exp(-Lc/Lp_guess)))
		error = (e2e-e2e_avg)/e2e_avg
		Lp_guess = (1-error)*Lp_guess
		counter += 1
		if counter > max_iter:
			print("Error: maximum number of iterations reached when calculating Lp from E2E distribution.")
			break

	### plotting
	plt.figure('E2E')
	plt.hist(e2e_arr,bins=e2e_bins)
	plt.axvline(e2e_avg,color="k")
	plt.xlabel('R')
	plt.ylabel('count')
	plt.title('True End-to-End Distance Distribution\nContour Length = ' + str(precision_round(Lc,2)) + " nm")
	plt.savefig('figures/e2e_dist')
	return Lp_guess


def theta_dist(v_tan):
	theta = np.zeros((int(v_tan.size/3),1))
	steps_use = int(v_tan.shape[0])
	np_scaf = int(v_tan.shape[1])
	for i in range(steps_use):
		for j in range(np_scaf):
			if j < np_scaf-2:
				u1 = v_tan[i,j,:]/np.linalg.norm(v_tan[i,j,:])
				u2 = v_tan[i,j+1,:]/np.linalg.norm(v_tan[i,j+1,:])
				dot = np.dot(u1,u2)
				if dot > 1:
					theta[i*np_scaf+j]
				else:
					theta[i*np_scaf+j] = math.acos(dot)
			else:
				theta[i*np_scaf+j] = theta[i*np_scaf+j-1]
	theta_avg = np.mean(theta)
	plt.figure('Theta')
	plt.hist(theta,bins=20)
	plt.axvline(theta_avg,color="k")
	plt.xlabel(r'$\theta$')
	plt.ylabel('count')
	plt.title('Theta Distribution')
	plt.savefig('figures/theta_dist')


def calc_MSD(points,dbox,n_bin,sim_time_ms):
	steps_use = int(points.shape[0])
	nstep_bin = int(np.floor(steps_use/n_bin))
	time = np.linspace(0,sim_time_ms/n_bin,nstep_bin+1)
	time = time[0:-1]
	dis2 = np.zeros((nstep_bin,n_bin))
	for q in range(n_bin):
		com0 = np.mean(points[q*nstep_bin,:,:],axis=0)
		com_prev = com0
		pbc = np.zeros(3)
		for i in range(1,nstep_bin):
			com = np.mean(points[q*nstep_bin+i,:,:],axis=0)
			for k in range(3):
				pbc[k] = pbc[k] + int(round((com[k]-com_prev[k])/dbox))
				dis2[i,q] += (com[k]-pbc[k]*dbox-com0[k])**2
			com_prev = com
	MSD = np.mean(dis2,axis=1)
	plt.figure('MSD')
	plt.plot(time,MSD)
	for q in range(min(n_bin,10)):
		plt.plot(time,dis2[:,q],'k',alpha=0.3)
	plt.xlabel('t [ms]')
	plt.ylabel('MSD [nm]')
	plt.title('Mean Square Displacement')
	plt.savefig('figures/msd')


def calc_persistence(datFile,nba_skip,coarse_strand,steps_skip,coarse_time,MSD_bins,e2e_bins):
	### read file
	np_scaf,steps_use,v_tan,points,dbox,sim_time_ms = read_dnafold(datFile,nba_skip,coarse_strand,steps_skip,coarse_time)

	### end-to-end parameters
	e2e_first = 0 			# chopped/coarsened bead index (start at 0) to start from in E2E calculations
	e2e_last = np_scaf-1	# chopped/coarsened bead index (start at 0) to end on in E2E calculations

	### calculate values along chain
	dot_avg_arr,ldot_sum_arr,e2e2_avg_arr,l0_avg = calc_dot_avg(v_tan,points,np_scaf,steps_use)
	lc_arr_nm = np.arange(np_scaf)*l0_avg

	### trajectory data analysis
	theta_dist(v_tan)
	lp_e2e = e2e_dist(points,l0_avg,e2e_first,e2e_last,e2e_bins)
	calc_MSD(points,dbox,MSD_bins,sim_time_ms)

	### estimate persistence length
	lp_corr = lp_from_correlation(lc_arr_nm,dot_avg_arr)
	lp_sum = ldot_sum_arr[-1]
	lp_bam = lc_arr_nm[1]/(1-dot_avg_arr[1])

	### print estimated persistence lengths
	print("Lp from tangent correlation fit (the literature standard:                    {:1.2f} nm".format(lp_corr))
	print("Lp from sum of tangent correlation (perfect for infinitely long strands):    {:1.2f} nm".format(lp_sum))
	print("Lp from mean bond angle (best for short strands):                            {:1.2f} nm".format(lp_bam))
	print("Lp from mean end-to-end distance (best for short simulations):               {:1.2f} nm".format(lp_e2e))

	### return calculations
	lp_est = lp_corr
	return lc_arr_nm,dot_avg_arr,ldot_sum_arr,e2e2_avg_arr,sim_time_ms,lp_est


def unit_vector(vector):
	return vector / np.linalg.norm(vector)


def calc_dot(v0,v1):
	v0_u = unit_vector(v0)
	v1_u = unit_vector(v1)
	dot = np.dot(v0_u, v1_u)
	if dot > 1:
		dot = 1
	return dot


def applyPBC(r,dbox):
	return r - dbox*np.round(r/dbox)


def calc_l0_avg(points):
	steps_use = int(points.shape[0])
	np_scaf = int(points.shape[1])
	max_l0 = np.zeros(steps_use)
	l0_sum = 0
	for i in range(steps_use):
		for j in range(np_scaf-1):
			l = np.linalg.norm(points[i,j+1,:] - points[i,j,:])
			l0_sum += l
			if l > max_l0[i]:
				max_l0[i] = l
		# print("At coarse step {}, max_l0 = {}".format(i,max_l0[i]))
	# print("Overall max_l0 = {}".format(np.max(max_l0)))
	l0_avg = l0_sum/(steps_use*(np_scaf-1))
	return l0_avg


def precision_round(number, digits):
	power = "{:e}".format(number).split('e')[1]
	return round(number,-(int(power)-digits))


### run the script
if __name__ == "__main__":
	main()