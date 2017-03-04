#!/bin/tcsh -f

if ( "$#" != "4" ) then
    echo "Missing parameters: printMDSTTree.sh mdstFileName runNo evtNo nEvts"
    exit -1
endif

source /sw/belle/local/etc/cshrc_general
setenv BASF_MODULE_DIR .:$BASF_MODULE_DIR

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
module put_parameter printtreemodule strVars\MXYZEi
module put_parameter printtreemodule runNo\\$run
module put_parameter printtreemodule evtNo\\$event
module put_parameter printtreemodule nEvts\\$n

initialize

process_event $mdst 0 

terminate

EOF
