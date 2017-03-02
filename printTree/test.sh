#!/bin/tcsh -f

if ( "$#" != "4" ) then
    echo "Missing parameters: printMDSTTree.sh mdstFileName runNo evtNo nEvts"
    exit -1
endif

setenv USE_GRAND_REPROCESS_DATA 1
source /sw/belle/local/etc/cshrc_general
source /home/belle/dcervenk/Env/cshrc.post
setenv BASF_MODULE_DIR ./bin:$BASF_MODULE_DIR

# Parameters
set mdst="$1"
set oFile="`basename $mdst .mdst`.txt"
set run="$2"
set event="$3"
set n="$4"

basf << EOF

path create main
path add_module main printtreemodule

path add_condition main >:0:EXIT
path add_condition main <=:0:KILL

module put_parameter printtreemodule fileName\\$oFile
module put_parameter printtreemodule qqOnly\1
module put_parameter printtreemodule strVars\MXYZE
module put_parameter printtreemodule runNo\\$run
module put_parameter printtreemodule evtNo\\$event
module put_parameter printtreemodule nEvts\\$n

initialize

process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/mixed/00/all/0127/on_resonance/09/evtgen-mixed-00-all-e000055r000937-b20090127_0910.mdst 0

terminate

EOF
