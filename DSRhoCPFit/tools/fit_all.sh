#!/usr/bin/zsh

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

DIR="cr_ti_Kpi_eff_190531"; mkdir results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config.json --fit=CR --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done

wait_till_all_processes_end DSRhoCPFit
echo "CR TI Kpi done"

DIR="crscf_ti_Kpi_eff_190531_scf_190531"; mkdir results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config.json --fit=CRSCF --scf-histo=scf_Kpi_190531.root --time-independent --log --cpus=1 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done

wait_till_all_processes_end DSRhoCPFit
echo "CRSCF TI Kpi done"

DIR="all_ti_Kpi_eff_190531_scf_190531"; mkdir results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config.json --fit=all --scf-histo=scf_Kpi_190531.root --time-independent --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).root ../data/Kpi/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done

wait_till_all_processes_end DSRhoCPFit
echo "All TI Kpi done"

DIR="cr_td_Kpi_eff_190531"; mkdir results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config.json --fit=CR --mixing --log --cpus=4 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done

wait_till_all_processes_end DSRhoCPFit
echo "CR TD Kpi done"

DIR="crscf_td_Kpi_eff_190531_scf_190531"; mkdir results/$DIR logs/$DIR; for FILE in ../data/DSRho-mdst_basf2_mod_real_unmod_<1-25>.root ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config.json --fit=CRSCF --scf-histo=scf_Kpi_190531.root --mixing --log --cpus=4 results/$DIR/$(basename -s.root $FILE).result $FILE &> logs/$DIR/$(basename -s.root $FILE).log &; done

wait_till_all_processes_end DSRhoCPFit
echo "CRSCF TD Kpi done"

DIR="all_td_Kpi_eff_190531_scf_190531"; mkdir results/$DIR logs/$DIR; for I in $(seq 0 5) ; do nice ./DSRhoCPFit --efficiency-model=6 --efficiency-file=eff_Kpi_190531.root --config=config.json --fit=all --scf-histo=scf_Kpi_190531.root --mixing --log --cpus=4 results/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).result ../data/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).root ../data/Kpi/mc_wo_signal/DSRhoSkim_svd*$I.root &> logs/$DIR/DSRho-mdst_basf2_mod_real_unmod_$((I+1)).log &; done

wait_till_all_processes_end DSRhoCPFit
echo "All TD Kpi done"

