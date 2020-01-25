#compdef DSRhoCPFit

# ----------------------------------------------------------------------
# zsh completions for dst-rho-ana programmes
#
# source this file in your .zshrc file to enable.
# ----------------------------------------------------------------------

_dsrhocpfit_complete()
{
	_arguments \
	'--cpus=[number of CPU cores to use for fitting and plotting]' \
	'--components=[fit the specified components]:type:->components' \
	'--config=[read in configuration from the specified file]:filename:_files' \
	'--events=[total number of events to be read from input files]' \
	'--exclude-channels=[exclude channels from fit]:type:->channels' \
	'--MC=[whether the fit is MC or data]:type:->binary' \
	'--version[show version]' \
	'--help[display help]' \
	'--time-independent[make a time-independent fit]' \
	'--log[save copy of log to results file]' \
	'--plot-dir=[create lifetime/mixing plots]:directory:_files' \
	'--perfect-tag[use MC info to get perfect tagging]' \
	'--fix=[fix specified argument(s) to input values in the fit]:parameter:->parameters' \
	'--generator-level[do a generator level fit]' \
	'--output=[output filename]:filename:_files'

	case "$state" in
		components)
			_values -s ' ' 'type' CR CRSCF all
			;;
		binary)
			_values -s ' ' 'type' 0 1
			;;
		channels)
			_values -s ' ' 'type' Kpi Kpipi0 K3pi
			;;
		jsonfiles)
			local -a json_files
			json_files=(*.json)
			_multi_parts / json_files
			;;
		rootfiles)
			local -a root_files
			root_files=(*.root)
			_multi_parts / root_files
			;;
		parameters)
			_values -s ',' 'parameters' all xy trans nota0 ap apa a0 ata xp x0 xt yp y0 yt xpb x0b xtb ypb y0b ytb
			;;
	esac
}

_dsrhobackground_complete()
{
	_arguments \
	'--cpus=[number of CPU cores to use for fitting and plotting]' \
	'--kde[use kernel density estimation]' \
	'--help[display help]' \
	'--physics[use a physics-based dt PDF]' \
	'--histo[create a histogram PDF]' \
	'--plot-dir=[create lifetime/mixing plots]:directory:_files' \
	'*:input files:_files'
}

_dsrholifetime_complete()
{
	_arguments \
	'--cpus=[number of CPU cores to use for fitting and plotting]' \
	'--components=[fit the specified components]:type:->components' \
	'--config=[read in configuration from the specified file]:filename:_files' \
	'--events=[total number of events to be read from input files]' \
	'--lifetime[do a lifetime fit]' \
	'--mixing[do a mixing fit]' \
	'--physics[use a physics-based BKG PDF]' \
	'--version[show version]' \
	'--help[display help]' \
	'--plot-dir=[create lifetime/mixing plots]:directory:_files' \
	'--perfect-tag[use MC info to get perfect tagging]'

	case "$state" in
		components)
			_values -s ' ' 'type' CR CRSCF all
			;;
	esac
}

_dsrhoyield_complete()
{
	_arguments \
	'--cpus=[number of CPU cores to use for fitting and plotting]' \
	'--MC[whether the fit is MC or data]' \
	'--version[show version]' \
	'--help[display help]' \
	'--plot-dir=[create lifetime/mixing plots]:directory:_files'
}

compdef _dsrhocpfit_complete DSRhoCPFit
compdef _dsrhobackground_complete DSRhoBackground
compdef _dsrholifetime_complete DSRhoLifetime
compdef _dsrhoyield_complete DSRhoYield

