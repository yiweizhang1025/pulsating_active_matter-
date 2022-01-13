gcc ep_inc_Param_order.c -o Order_ep_inc -lm
#"Input: output Lx Ly rho T mu_r mu_phi omega dt totaltime initial_time interval_record interval_snap seed radius amplitude_r amplitude_phi_max amplitude_phi_min lambda dist drho dx dphi"
#align_list=('62' '64' '66' '68' '70' '72' '74' '76' '78' '80' '82' '84' '86' '88' '90')
#align_list=('1') # '30' '50' '80' '82' '84' '86' '88' '90')
traj_list=('2' '3' '4' '5' '6' '7' '8' '9' '10') # '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37' '38' '39' '40' '41' '42' '43' '44' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54' '55' '56' '57' '58' '59' '60' '61' '62' '63' '64' '65' '66' '67' '68' '69' '70' '71' '72' '73' '74' '75' '76' '77' '78' '79' '80' '81' '82' '83' '84' '85' '86' '87' '88' '89' '90' '91' '92' '93' '94' '95' '96' '97' '98' '99' '100')
export nb_traj=${#traj_list[@]}

for dens in '0.8'
do
#	export dens
#	for align in ${align_list[@]}
#       	do
	 ### Check if a directory does not exist ###
 #       if [ ! -d "align_$align" ]
 #       then
  #             mkdir align_"$align"
 #       fi
 #       rm  align_"$align"/*
		for traj in "${traj_list[@]}"
		do
        	./Order_ep_inc ep_inc_adp_rho_"$dens"_omega_1e1_amp_1e0_alignmax_1.5e1_alignmin_0e1_lambda_5e-2_trj_"$traj" 1e2 1e2 "$dens"e0 1e0 1e0 1e0 1e1 5e-3 160e2 12e2 1e-1 1e-1 $RANDOM 1e0 1e0 1.5e1 0e1 5e-2 8e0 1e-3 1e-2 1e-2 &
		done 
	wait
	done
done
python writelines.py $nb_traj
source ./avg.sh
gnuplot plot.pl
