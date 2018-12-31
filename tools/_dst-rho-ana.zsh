#compdef DSRhoCPFit

# ----------------------------------------------------------------------
# zsh completions for dst-rho-ana programmes
#
# source this file in your .zshrc file to enable.
# ----------------------------------------------------------------------

_dsrhocpfit_complete()
{
	_arguments '--cpus=[number of CPU cores to use for fitting and plotting]' '--efficiency-model=[number of the efficiency model to be used]:model:->effmodel' '--fit=[do a specified fit type]:type:->fittype' '--config=[read in configuration from the specified file]:filename:->jsonfiles''--help[display help]' '--time-independent[make a time-independent fit]' '--log[save copy of log to results file]' '--mixing[make a mixing fit]' '--events=[number of events to be imported from the input file]' '--plot-dir=[create lifetime/mixing plots]:directory:_files' '--perfect-tag[use MC info to get perfect tagging]' '--fix=[fix specified argument(s) to input values in the fit]:parameter:->parameters'
	case "$state" in
		effmodel)
			_values -s ' ' 'model' 0 1 2 3 4 5 6
			;;
		fittype)
			_values -s ' ' 'type' CR CRSCF all
			;;
		jsonfiles)
			local -a json_files
			json_files=(*.json)
			_multi_parts / json_files
			;;
		parameters)
			_values -s ',' 'parameters' all xy trans nota0 ap apa a0 ata xp x0 xt yp y0 yt xpb x0b xtb ypb y0b ytb
			;;
	esac
}

compdef _dsrhocpfit_complete DSRhoCPFit
