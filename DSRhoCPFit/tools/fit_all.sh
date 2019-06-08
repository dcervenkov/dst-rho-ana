#!/usr/bin/zsh
# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o nounset

# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null 
if [[ ${pipestatus[1]} -ne 4 ]]; then
	echo 'I’m sorry, `getopt --test` failed in this environment.'
	exit 1
fi

OPTIONS=pntdicsa
LONGOPTS=kpi,kpipi0,k3pi,td,ti,cr,crscf,all

# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${pipestatus[1]} -ne 0 ]]; then
	# e.g. return value is 1
	#  then getopt has complained about wrong arguments to stdout
	exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

KPI=0 KPIPI0=0 K3PI=0 TD=0 TI=0 CR=0 CRSCF=0 ALL=0
# now enjoy the options in order and nicely split until we see --
while true; do
	case "$1" in
		-p|--kpi)
			KPI=1
			shift
			;;
		-n|--kpipi0)
			KPIPI0=1
			shift
			;;
		-t|--k3pi)
			K3PI=1
			shift
			;;
		-d|--td)
			TD=1
			shift
			;;
		-i|--ti)
			TI=1
			shift
			;;
		-c|--cr)
			CR=1
			shift
			;;
		-s|--crscf)
			CRSCF=1
			shift
			;;
		-a|--all)
			ALL=1
			shift
			;;
		--)
			shift
			break
			;;
		*)
			echo "Unknown argument"
			exit 3
			;;
	esac
done

wait_till_all_processes_end() {
PROC_NAME=$1
while true; do
	if pgrep -x $PROC_NAME > /dev/null; then
		FOUND=true
		sleep 10
	else
		return 0
	fi
done
}


cd ..

if [ "$TI" = 1 ]; then
	if [ "$KPI" = 1 ]; then
		if [ "$CR" = 1 ]; then
			DIR="cr_ti_Kpi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config_Kpi.json --fit=CR --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CR TI Kpi done"
		fi

		if [ "$CRSCF" = 1 ]; then
			DIR="crscf_ti_Kpi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config_Kpi.json --fit=CRSCF --scf-histo=scf_Kpi_190531.root --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CRSCF TI Kpi done"
		fi

		if [ "$ALL" = 1 ]; then
			DIR="all_ti_Kpi"; mkdir -p results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config_Kpi.json --fit=all --scf-histo=scf_Kpi_190531.root --time-independent --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).root ../data/Kpi/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "All TI Kpi done"
		fi
	fi

	if [ "$KPIPI0" = 1 ]; then
		if [ "$CR" = 1 ]; then
			DIR="cr_ti_Kpipi0"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190528_Kpipi0/DSRho-mdst_basf2_190528_<001-100>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpipi0_190603.root --config=config_Kpipi0.json --fit=CR --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CR TI Kpipi0 done"
		fi

		if [ "$CRSCF" = 1 ]; then
			DIR="crscf_ti_Kpipi0"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190528_Kpipi0/DSRho-mdst_basf2_190528_<001-100>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpipi0_190603.root --config=config_Kpipi0.json --fit=CRSCF --scf-histo=scf_Kpipi0_190603.root --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CRSCF TI Kpipi0 done"
		fi

		if [ "$ALL" = 1 ]; then
			DIR="all_ti_Kpipi0"; mkdir -p results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpipi0_190603.root --config=config_Kpipi0.json --fit=all --scf-histo=scf_Kpipi0_190603.root --time-independent --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/basf2_190528_Kpipi0/DSRho-mdst_basf2_190528_00$((I+1))_svd2.root ../data/Kpipi0/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "All TI Kpipi0 done"
		fi
	fi

	if [ "$K3PI" = 1 ]; then
		if [ "$CR" = 1 ]; then
			DIR="cr_ti_K3pi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190529_K3pi/DSRho-mdst_basf2_190529_<001-100>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_K3pi_190603.root --config=config_K3pi.json --fit=CR --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CR TI K3pi done"
		fi

		if [ "$CRSCF" = 1 ]; then
			DIR="crscf_ti_K3pi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190529_K3pi/DSRho-mdst_basf2_190529_<001-100>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_K3pi_190603.root --config=config_K3pi.json --fit=CRSCF --scf-histo=scf_K3pi_190603.root --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CRSCF TI K3pi done"
		fi

		if [ "$ALL" = 1 ]; then
			DIR="all_ti_K3pi"; mkdir -p results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_K3pi_190603.root --config=config_K3pi.json --fit=all --scf-histo=scf_K3pi_190603.root --time-independent --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/basf2_190529_K3pi/DSRho-mdst_basf2_190529_00$((I+1))_svd2.root ../data/K3pi/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "All TI K3pi done"
		fi
	fi
fi

if [ "$TD" = 1 ]; then
	if [ "$KPI" = 1 ]; then
		if [ "$CR" = 1 ]; then
			DIR="cr_td_Kpi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config_Kpi.json --fit=CR --mixing --log --cpus=4 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CR TD Kpi done"
		fi

		if [ "$CRSCF" = 1 ]; then
			DIR="crscf_td_Kpi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config_Kpi.json --fit=CRSCF --scf-histo=scf_Kpi_190531.root --mixing --log --cpus=4 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CRSCF TD Kpi done"
		fi

		if [ "$ALL" = 1 ]; then
			DIR="all_td_Kpi"; mkdir -p results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config_Kpi.json --fit=all --scf-histo=scf_Kpi_190531.root --mixing --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/DSRho-mdst_basf2_mod_real_unmod_00$((I+1)).root ../data/Kpi/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "All TD Kpi done"
		fi

	fi

	if [ "$KPIPI0" = 1 ]; then
		if [ "$CR" = 1 ]; then
			DIR="cr_td_Kpipi0"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190528_Kpipi0/DSRho-mdst_basf2_190528_<001-025>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpipi0_190603.root --config=config_Kpipi0.json --fit=CR --mixing --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CR TD Kpipi0 done"
		fi

		if [ "$CRSCF" = 1 ]; then
			DIR="crscf_td_Kpipi0"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190528_Kpipi0/DSRho-mdst_basf2_190528_<001-025>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpipi0_190603.root --config=config_Kpipi0.json --fit=CRSCF --scf-histo=scf_Kpipi0_190603.root --mixing --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CRSCF TD Kpipi0 done"
		fi

		if [ "$ALL" = 1 ]; then
			DIR="all_td_Kpipi0"; mkdir -p results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpipi0_190603.root --config=config_Kpipi0.json --fit=all --scf-histo=scf_Kpipi0_190603.root --mixing --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/basf2_190528_Kpipi0/DSRho-mdst_basf2_190528_00$((I+1))_svd2.root ../data/Kpipi0/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "All TD Kpipi0 done"
		fi
	fi

	if [ "$K3PI" = 1 ]; then
		if [ "$CR" = 1 ]; then
			DIR="cr_td_K3pi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190529_K3pi/DSRho-mdst_basf2_190529_<001-025>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_K3pi_190603.root --config=config_K3pi.json --fit=CR --mixing --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CR TD K3pi done"
		fi

		if [ "$CRSCF" = 1 ]; then
			DIR="crscf_td_K3pi"; mkdir -p results/$DIR logs/$DIR; for FILE in ../data/basf2_190529_K3pi/DSRho-mdst_basf2_190529_<001-025>_svd2.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_K3pi_190603.root --config=config_K3pi.json --fit=CRSCF --scf-histo=scf_K3pi_190603.root --mixing --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "CRSCF TD K3pi done"
		fi

		if [ "$ALL" = 1 ]; then
			DIR="all_td_K3pi"; mkdir -p results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_K3pi_190603.root --config=config_K3pi.json --fit=all --scf-histo=scf_K3pi_190603.root --mixing --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/basf2_190529_K3pi/DSRho-mdst_basf2_190529_$((I+1))_svd2.root ../data/K3pi/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done
			wait_till_all_processes_end DSRhoCPFit
			echo "All TD K3pi done"
		fi
	fi
fi
