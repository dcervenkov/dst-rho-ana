import glob
import shutil

STREAMS = range(0, 6)
RBINS = range(0, 7)
CHANNELS = ["Kpi", "Kpipi0", "K3pi"]
CHANNELS_AND_TOGETHER = ["Kpi", "Kpipi0", "K3pi", "together"]
SVDS = [1, 2]
BKG_TYPES = ["scf", "bkg"]

MC_SIGNAL_6_STREAMS = {
    "Kpi": glob.glob("data/Kpi/signal_mc/*_basf2_0[0-5]_*.root"),
    "Kpipi0": glob.glob("data/Kpipi0/signal_mc/*_basf2_0[0-5]_*.root"),
    "K3pi": glob.glob("data/K3pi/signal_mc/*_basf2_0[0-5]_*.root")
}
MC_WO_SIGNAL = {
    "Kpi": glob.glob("data/Kpi/mc_wo_signal/*.root"),
    "Kpipi0": glob.glob("data/Kpipi0/mc_wo_signal/*.root"),
    "K3pi": glob.glob("data/K3pi/mc_wo_signal/*.root")
}
DATA_SIDEBANDS = {
    "Kpi": glob.glob("data/Kpi/sidebands/*.root"),
    "Kpipi0": glob.glob("data/Kpipi0/sidebands/*.root"),
    "K3pi": glob.glob("data/K3pi/sidebands/*.root")
}

wildcard_constraints:
    randomization = "[a-zA-Z0-9_]*"


def get_data_for_background(wildcards):
    channel = wildcards.channel
    mc = wildcards.mc
    bkg_type = wildcards.bkg_type
    if channel == "together":
        if mc == "mc":
            if bkg_type == "scf":
                return (MC_SIGNAL_6_STREAMS["Kpi"] +
                        MC_SIGNAL_6_STREAMS["Kpipi0"] +
                        MC_SIGNAL_6_STREAMS["K3pi"])
            elif bkg_type == "bkg":
                return (MC_WO_SIGNAL["Kpi"] +
                        MC_WO_SIGNAL["Kpipi0"] +
                        MC_WO_SIGNAL["K3pi"])
            elif bkg_type == "mixed":
                return glob.glob("data/K*/mc_wo_signal/*mixed*.root")
            elif bkg_type == "charged":
                return glob.glob("data/K*/mc_wo_signal/*charged*.root")
            elif bkg_type == "charmuds":
                files = glob.glob("data/K*/mc_wo_signal/*charm*.root")
                files.extend(glob.glob("data/K*/mc_wo_signal/*uds*.root"))
                return files

        else:
            return (DATA_SIDEBANDS["Kpi"] +
                    DATA_SIDEBANDS["Kpipi0"] +
                    DATA_SIDEBANDS["K3pi"])
    else:
        if mc == "mc":
            if bkg_type == "scf":
                return MC_SIGNAL_6_STREAMS[channel]
            else:
                return MC_WO_SIGNAL[channel]
        else:
            return DATA_SIDEBANDS[channel]


def get_params_for_background(wildcards):
    bkg_type = wildcards.bkg_type
    if bkg_type == "scf":
        return "--notail --nodelta --physics"
    elif bkg_type == "bkg":
        return "--notail --physics"
    elif bkg_type == "sidebands":
        return "--notail --physics"
    elif bkg_type == "mixed":
        return "--notail --wtag --mixing"
    elif bkg_type == "charged":
        return "--notail --wtag"
    elif bkg_type == "charmuds":
        return "--physics"


def get_excluded_channels(wildcards):
    if wildcards.channel == "together":
        return ""
    else:
        excluded_channels = set(CHANNELS) - set(wildcards.channel.split("-"))
        if excluded_channels:
            return "--exclude-channels=" + ",".join(excluded_channels)
        else:
            return ""

def create_stream_config(old_stream, old_file, new_stream, new_file):
    new_config = ""
    with open(old_file, "r") as f:
        for line in f:
            new_config += line.replace(
                f"stream{old_stream}", f"stream{new_stream}"
            ).replace(f"_{int(old_stream):02}_", f"_{int(new_stream):02}_"
                    ).replace(f"{old_stream}.root", f"{new_stream}.root")
    with open(new_file, "w") as f:
        f.write(new_config)

def create_rnd_config(old_file, new_file, number, line_tag, prefix):
    new_config = ""
    with open(old_file, "r") as f:
        for line in f:
            if line_tag in line:
                new_config += line.replace(prefix, "randomized/" + prefix
                    ).replace(".root", f"_rnd_{number}.root")
            else:
                new_config += line
    with open(new_file, "w") as f:
        f.write(new_config)

rule all:
    input:
        yield_jobs = (
            expand("DSRhoYield/plots/{channel}_stream{stream}/fit_results.root",
                   channel=CHANNELS, stream=STREAMS),
            expand("DSRhoYield/plots/{channel}/fit_results.root", channel=CHANNELS),

            expand("DSRhoYield/plots/{channel}_stream{stream}_rbin{rbin}/fit_results.root",
                   channel=CHANNELS, stream=STREAMS, rbin=RBINS),
            expand("DSRhoYield/plots/{channel}_rbin{rbin}/fit_results.root", channel=CHANNELS, rbin=RBINS)
            ),

        background_jobs = (
            expand("DSRhoBackground/results/{channel}_mc_{bkg_type}.json",
                   channel=["together"], bkg_type=["mixed", "charged", "charmuds"]),

            expand("DSRhoBackground/results/{channel}_mc_{bkg_type}.json",
                   channel=CHANNELS_AND_TOGETHER, bkg_type=BKG_TYPES),
            expand("DSRhoBackground/results/{channel}_data_sidebands.json",
                   channel=CHANNELS_AND_TOGETHER),

            expand("DSRhoBackground/results/nonphys_{channel}_mc_{bkg_type}.json",
                   channel=CHANNELS_AND_TOGETHER, bkg_type=BKG_TYPES),
            expand("DSRhoBackground/results/nonphys_{channel}_data_sidebands.json",
                   channel=CHANNELS_AND_TOGETHER)
            ),
        lifetime_jobs = (
            expand("DSRhoLifetime/results/{channel}/{component}_{type}_{stream}",
                   channel=CHANNELS_AND_TOGETHER, component=["CR", "CRSCF"], type=["lifetime", "mixing"],
                   stream=[f"{stream:02}" for stream in range(99)]),
            expand("DSRhoLifetime/results/{channel}/{component}_{type}_{stream}",
                   channel=CHANNELS_AND_TOGETHER, component=["all"], type=["lifetime", "mixing"],
                   stream=range(6)),

            expand("DSRhoLifetime/results/{channel}_{component}_{type}_{stream}",
                   channel=CHANNELS_AND_TOGETHER, component=["all"], type=["lifetime", "mixing"],
                   stream=range(6)),
            expand("DSRhoLifetime/plots/{channel}_{component}_{type}_{stream}",
                   channel=CHANNELS_AND_TOGETHER, component=["all"], type=["lifetime", "mixing"],
                   stream=range(6)),

            expand("DSRhoLifetime/results/{channel}_{type}_data",
                   channel=CHANNELS_AND_TOGETHER, component=["all"], type=["lifetime", "mixing"]),
            expand("DSRhoLifetime/plots/{channel}_{type}_data",
                   channel=CHANNELS_AND_TOGETHER, component=["all"], type=["lifetime", "mixing"])
            ),
        cpfit_jobs = (
            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
            #     channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["CR", "CRSCF"],
            #     stream=range(99), group=["nominal"]),
            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
            #     channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["all"],
            #     stream=range(6), group=["nominal"]),

            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}_plot/stream{stream}",
            #     channel=CHANNELS_AND_TOGETHER, type=["ti"], components=["CR", "CRSCF"],
            #     stream=range(1), group=["nominal"]),
            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}_plot/stream{stream}",
            #     channel=CHANNELS_AND_TOGETHER, type=["ti"], components=["all"],
            #     stream=range(1), group=["nominal"]),

            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}",
            #     channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg", "sidebkg"],
            #     group=["nominal"]),
            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}",
            #     channel=["Kpi-K3pi"], type=["ti", "td"], bkg=["mcbkg", "sidebkg"],
            #     group=["nominal"])

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["CR", "CRSCF"],
                stream=range(99), group=["nominal_rbin"]),
            expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["all"],
                stream=range(6), group=["nominal_rbin"]),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}_plot/stream{stream}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["CR", "CRSCF"],
                stream=range(1), group=["nominal_rbin"]),
            expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}_plot/stream{stream}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["all"],
                stream=range(1), group=["nominal_rbin"]),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg", "sidebkg"],
                group=["nominal_rbin"]),
            ),
        cpfit_extra_jobs = (
            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
            #     channel=["Kpi"], type=["ti", "td"], components=["CR", "CRSCF"], stream=range(8),
            #     group=["correct_phi3"]),
            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
            #     channel=["Kpi"], type=["ti", "td"], components=["all"], stream=range(6),
            #     group=["correct_phi3"]),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_eff"], configno=range(100)),
            expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["all"], stream=range(6),
                group=["nonres"]),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["all"], stream=range(6),
                group=["nonres_mcflag7"]),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], components=["all"], stream=range(6),
                group=["nonres_0.3_mcflag7"]),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_scf"], configno=range(99)),

            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
            #     channel=CHANNELS_AND_TOGETHER, type=["td"], bkg=["mcbkg"],
            #     group=["randomized_scf_dt"], configno=range(99)),

            # expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
            #     channel=CHANNELS_AND_TOGETHER, type=["td"], bkg=["mcbkg"],
            #     group=["randomized_scf_dt_corr"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_bkg"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_bkg_corr"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["td"], bkg=["mcbkg"],
                group=["randomized_bkg_dt"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["td"], bkg=["mcbkg"],
                group=["randomized_bkg_dt_corr"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_yield"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_yield_corr"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_all"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["ti", "td"], bkg=["mcbkg"],
                group=["randomized_all_corr"], configno=range(99)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=CHANNELS_AND_TOGETHER, type=["td"], bkg=["mcbkg"],
                group=["sigma_tau_dm"], configno=range(5)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=["Kpipi0"], type=["ti", "td"], bkg=["mcbkg"],
                group=["measured_scf"], configno=range(1)),

            expand("DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
                channel=["Kpi", "Kpipi0"], type=["ti", "td"], bkg=["mcbkg"],
                group=["scf_swap"], configno=range(1)),
            )

rule yield_jobs:
    input:
        rules.all.input.yield_jobs

rule background_jobs:
    input:
        rules.all.input.background_jobs

rule lifetime_jobs:
    input:
        rules.all.input.lifetime_jobs

rule cpfit_jobs:
    input:
        rules.all.input.cpfit_jobs,
        rules.all.input.cpfit_extra_jobs

rule yield_mc:
    input:
        "data/{channel}/realistic_mc",
        "data/{channel}/realistic_mc/stream{stream}"
    output:
        "DSRhoYield/plots/{channel}_stream{stream}/fit_results.root"
    log:
        "DSRhoYield/logs/{channel}_stream{stream}"
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together",
        stream = "\d+"
    run:
        shell("./DSRhoYield/DSRhoYield --MC {input} " + os.path.dirname(output[0]) + " &> {log}")

rule yield_data:
    input:
        "data/{channel}/realistic_mc",
        "data/{channel}"
    output:
        "DSRhoYield/plots/{channel}/fit_results.root"
    log:
        "DSRhoYield/logs/{channel}"
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together"
    run:
        shell("./DSRhoYield/DSRhoYield {input} " + os.path.dirname(output[0]) + " &> {log}")

rule yield_mc_rbin:
    input:
        "data/{channel}/realistic_mc",
        "data/{channel}/realistic_mc/stream{stream}"
    output:
        "DSRhoYield/plots/{channel}_stream{stream}_rbin{rbin}/fit_results.root"
    log:
        "DSRhoYield/logs/{channel}_stream{stream}_rbin{rbin}"
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together",
        stream = "\d+",
        rbin = "\d+"
    run:
        shell("./DSRhoYield/DSRhoYield --MC --rbin {wildcards.rbin} {input} " + os.path.dirname(output[0]) + " &> {log}")

rule yield_data_rbin:
    input:
        "data/{channel}/realistic_mc",
        "data/{channel}"
    output:
        "DSRhoYield/plots/{channel}_rbin{rbin}/fit_results.root"
    log:
        "DSRhoYield/logs/{channel}_rbin{rbin}"
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together",
        rbin = "\d+"
    run:
        shell("./DSRhoYield/DSRhoYield --rbin {wildcards.rbin} {input} " + os.path.dirname(output[0]) + " &> {log}")

rule yield_summary:
    input:
        expand("DSRhoYield/plots/{channel}_stream{stream}/fit_results.root",
               channel=CHANNELS, stream=STREAMS),
        expand("DSRhoYield/plots/{channel}/fit_results.root",
               channel=CHANNELS)
    output:
        expand("DSRhoYield/results/{channel}_{type}_fractions.json", channel=CHANNELS, type=["mc", "data"]),
        "DSRhoYield/results/avg_data_fractions.json"
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together",
        stream = "\d+",
    shell:
        "cd DSRhoYield && ./tools/print_all_yields.py &> yield.log"

rule yield_summary_rbin:
    input:
        expand("DSRhoYield/plots/{channel}_stream{stream}_rbin{{rbin}}/fit_results.root",
               channel=CHANNELS, stream=STREAMS),
        expand("DSRhoYield/plots/{channel}_rbin{{rbin}}/fit_results.root",
               channel=CHANNELS)
    output:
        expand("DSRhoYield/results/rbin{{rbin}}/{channel}_{type}_fractions.json",
               channel=CHANNELS, type=["mc", "data"]),
        "DSRhoYield/results/rbin{rbin}/avg_data_fractions.json"
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together",
        stream = "\d+",
        rbin = "\d+"
    shell:
        "cd DSRhoYield && ./tools/print_all_yields.py --directory results/rbin{wildcards.rbin} --rbin {wildcards.rbin} &> yield_rbin{wildcards.rbin}.log"

rule yield_summary_rbin_randomized:
    input:
        expand("DSRhoYield/plots/{channel}_rbin{{rbin}}/fit_results.root",
               channel=CHANNELS)
    output:
        expand("DSRhoYield/results_{{randomization}}/rbin{{rbin}}/{channel}_{type}_fractions.json",
               channel=CHANNELS, type=["data"]),
        "DSRhoYield/results_{randomization}/rbin{rbin}/avg_data_fractions.json",
        expand("DSRhoYield/results_{{randomization}}/rbin{{rbin}}/rnd_{rnd}/{channel}_{type}_fractions.json",
               channel=CHANNELS, type=["data"], rnd=range(100)),
        expand("DSRhoYield/results_{{randomization}}/rbin{{rbin}}/rnd_{rnd}/avg_data_fractions.json",
               rnd=range(100))
    log:
        "logs/yield_{randomization}_rbin{rbin}.log"
    params:
        "--randomize 100",
        lambda wildcards:
            "--correlations" if "corr" in wildcards.randomization else "",
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together",
        rbin = "\d+"
    shell:
        "cd DSRhoYield && ./tools/print_all_yields.py --directory results_{wildcards.randomization}/rbin{wildcards.rbin} {params} --rbin {wildcards.rbin} &> {log}"

rule yield_rearrange_randomized:
    input:
        expand("DSRhoYield/results_{{randomization}}/rbin{rbin}/rnd_{rnd}/{channel}_{type}_fractions.json",
               channel=CHANNELS, type=["data"], rbin=range(7), rnd=range(100)),
        expand("DSRhoYield/results_{{randomization}}/rbin{rbin}/rnd_{rnd}/avg_data_fractions.json",
               rbin=range(7), rnd=range(100)),
        expand("DSRhoYield/results/{channel}_{type}_fractions.json", channel=CHANNELS, type=["data"]),
        "DSRhoYield/results/avg_data_fractions.json"
    output:
        expand("DSRhoYield/results_{{randomization}}/rnd_{rnd}/rbin{rbin}/{channel}_{type}_fractions.json",
               channel=CHANNELS, type=["data"], rbin=range(7), rnd=range(100)),
        expand("DSRhoYield/results_{{randomization}}/rnd_{rnd}/rbin{rbin}/avg_data_fractions.json",
               rbin=range(7), rnd=range(100))
    run:
        for channel in [*CHANNELS, "avg"]:
            for rnd in range(100):
                for rbin in range(7):
                    shutil.move(f"DSRhoYield/results_{wildcards.randomization}/rbin{rbin}/rnd_{rnd}/{channel}_data_fractions.json",
                        f"DSRhoYield/results_{wildcards.randomization}/rnd_{rnd}/rbin{rbin}")
                shutil.copy(f"DSRhoYield/results/{channel}_data_fractions.json",
                    f"DSRhoYield/results_{wildcards.randomization}/rnd_{rnd}")
        for rbin in range(7):
            shutil.rmtree(f"DSRhoYield/results_{wildcards.randomization}/rbin{rbin}")

rule background:
    input:
        get_data_for_background
    output:
        result = "DSRhoBackground/results/{channel}_{mc}_{bkg_type}.json",
        plotdir = directory("DSRhoBackground/plots/{channel}_{mc}_{bkg_type}")
    log:
        "DSRhoBackground/logs/{channel}_{mc}_{bkg_type}"
    params:
        get_params_for_background
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together",
        mc = "mc|data"
    shell:
        "./DSRhoBackground/DSRhoBackground"
        " {params} --angular --plot-dir={output.plotdir} {output.result} {input} &> {log}"

rule background_nonphys:
    input:
        get_data_for_background
    output:
        result = "DSRhoBackground/results/nonphys_{channel}_{mc}_{bkg_type}.json",
        plotdir = directory("DSRhoBackground/plots/nonphys_{channel}_{mc}_{bkg_type}")
    log:
        "DSRhoBackground/logs/nonphys_{channel}_{mc}_{bkg_type}"
    shell:
        "./DSRhoBackground/DSRhoBackground"
        " --plot-dir={output.plotdir} --angular --empirical {output.result} {input} &> {log}"

rule lifetime_configs:
    input:
        expand(rules.background.output, channel=CHANNELS_AND_TOGETHER, mc="mc", bkg_type=BKG_TYPES),
        expand(rules.background.output, channel=CHANNELS_AND_TOGETHER, mc="data", bkg_type="sidebands"),
        rules.yield_summary.output,
        template = "DSRhoLifetime/configs/{group}/templates/{config}.template.json"
    output:
        "DSRhoLifetime/configs/{group}/{config}.json"
    wildcard_constraints:
        group = "[a-zA-Z0-9_]*"
    shell:
        "./tools/config_from_template.py {input.template} > {output}"

rule lifetime_no_plots_mc:
    input:
        config = lambda wildcards:
            "DSRhoLifetime/configs/config_mc_all.json" if wildcards.channel == "together" else
            "DSRhoLifetime/configs/config_mc.json",
        data = lambda wildcards:
            glob.glob(f"data/{wildcards.channel}/realistic_mc/stream{wildcards.stream}/*.root") if wildcards.component == "all" else
            glob.glob(f"data/{wildcards.channel}/signal_mc/*_{wildcards.stream}_svd?.root")

    output:
        "DSRhoLifetime/results/{channel}/{component}_{type}_{stream}"
    log:
        "DSRhoLifetime/logs/{channel}/{component}_{type}_{stream}"
    params:
        "--cpus=1 ",
        "--components={component} ",
        "--channel={channel}",
        "--physics",
        "--{type}"
    shell:
        "./DSRhoLifetime/DSRhoLifetime {params} --config={input.config} {output} {input.data} &> {log}"

rule lifetime_plots_mc:
    input:
        config = lambda wildcards:
            "DSRhoLifetime/configs/config_mc_all.json" if wildcards.channel == "together" else
            "DSRhoLifetime/configs/config_mc.json",
        data = lambda wildcards:
            glob.glob(f"data/{wildcards.channel}/realistic_mc/stream{wildcards.stream}/*.root") if wildcards.component == "all" else
            glob.glob(f"data/{wildcards.channel}/signal_mc/*_{wildcards.stream}_svd?.root")

    output:
        plotdir = directory("DSRhoLifetime/plots/{channel}_{component}_{type}_{stream}"),
        result = "DSRhoLifetime/results/{channel}_{component}_{type}_{stream}"
    log:
        "DSRhoLifetime/logs/{channel}_{component}_{type}_{stream}"
    params:
        "--cpus=1 ",
        "--components={component} ",
        "--channel={channel}",
        "--physics",
        "--{type}"
    shell:
        "./DSRhoLifetime/DSRhoLifetime {params} --config={input.config} "
        "--plot-dir={output.plotdir} "
        "{output.result} {input.data} &> {log}"

rule lifetime_plots_data:
    input:
        config = lambda wildcards:
            "DSRhoLifetime/configs/config_data_all.json" if wildcards.channel == "together" else
            "DSRhoLifetime/configs/config_data.json",
        data = lambda wildcards:
            glob.glob(f"data/{wildcards.channel}/*.root")
    output:
        plotdir = directory("DSRhoLifetime/plots/{channel}_{type}_data"),
        result = "DSRhoLifetime/results/{channel}_{type}_data"
    log:
        "DSRhoLifetime/logs/{channel}_{type}_data"
    params:
        "--cpus=1 ",
        "--components=all ",
        "--channel={channel}",
        "--physics",
        "--{type}"
    shell:
        "./DSRhoLifetime/DSRhoLifetime {params} --config={input.config} "
        "--plot-dir={output.plotdir} "
        "{output.result} {input.data} &> {log}"

rule cpfit_configs:
    input:
        expand(rules.background.output, channel=CHANNELS_AND_TOGETHER, mc="mc", bkg_type=BKG_TYPES),
        expand(rules.background_nonphys.output, channel=CHANNELS_AND_TOGETHER, mc="mc", bkg_type=BKG_TYPES),
        expand(rules.background.output, channel=CHANNELS_AND_TOGETHER, mc="data", bkg_type="sidebands"),
        expand(rules.background_nonphys.output, channel=CHANNELS_AND_TOGETHER, mc="data", bkg_type="sidebands"),
        rules.yield_summary.output,
        expand(rules.yield_summary_rbin.output, rbin=RBINS),
        template = "DSRhoCPFit/configs/{group}/templates/{config}.template.json"
    output:
        "DSRhoCPFit/configs/{group}/{config}.json"
    wildcard_constraints:
        group = "[a-zA-Z0-9_]*",
        config = "[a-zA-Z0-9_]*"
    shell:
        "./tools/config_from_template.py {input.template} > {output}"

rule cpfit_stream_configs:
    input:
        "DSRhoCPFit/configs/{group}/config_mc.json"
    output:
        "DSRhoCPFit/configs/{group}/streams/config_mc_{stream}.json"
    wildcard_constraints:
        stream = "\d+"
    run:
        create_stream_config(0, input[0], int(wildcards.stream), output[0])

rule cpfit_mc:
    input:
        config = "DSRhoCPFit/configs/{group}/streams/config_mc_{stream}.json",
    output:
        result = "DSRhoCPFit/results/{group}/{channel}_{type}_{components}/stream{stream}"
    log:
        "DSRhoCPFit/logs/{group}/{channel}_{type}_{components}/stream{stream}.log"
    wildcard_constraints:
        components = "CR|CRSCF|all"
    params:
        "--MC=1",
        "--cpus=1",
        "--log",
        "--components={components}",
        lambda wildcards:
            "--time-independent" if wildcards.type == "ti" else "",
        lambda wildcards:
            get_excluded_channels(wildcards)
    shell:
        "./DSRhoCPFit/DSRhoCPFit {params} --config={input.config} --output={output.result} &> {log}"

rule cpfit_mc_plot:
    input:
        config = "DSRhoCPFit/configs/{group}/streams/config_mc_{stream}.json",
    output:
        result = "DSRhoCPFit/results/{group}/{channel}_{type}_{components}_plot/stream{stream}",
        plotdir = directory("DSRhoCPFit/plots/{group}/{channel}_{type}_{components}/stream{stream}")
    log:
        "DSRhoCPFit/logs/{group}/{channel}_{type}_{components}_plot/stream{stream}.log"
    wildcard_constraints:
        components = "CR|CRSCF|all"
    params:
        "--MC=1",
        "--cpus=1",
        "--log",
        "--components={components}",
        lambda wildcards:
            "--time-independent" if wildcards.type == "ti" else "",
        lambda wildcards:
            get_excluded_channels(wildcards)
    shell:
        "./DSRhoCPFit/DSRhoCPFit {params} --config={input.config} --plot-dir={output.plotdir} "
        "--output={output.result} &> {log}"

rule cpfit_data_plot:
    input:
        config = "DSRhoCPFit/configs/{group}/config_data_{bkg}.json",
    output:
        result = "DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}",
        plotdir = directory("DSRhoCPFit/plots/{group}/{channel}_{type}_data_{bkg}")
    log:
        "DSRhoCPFit/logs/{group}/{channel}_{type}_data_{bkg}.log"
    wildcard_constraints:
        bkg = "[a-zA-z0-9_]*"
    params:
        "--MC=0",
        "--cpus=1",
        "--log",
        "--components=all",
        lambda wildcards:
            "--time-independent" if wildcards.type == "ti" else "",
        lambda wildcards:
            get_excluded_channels(wildcards)
    shell:
        "./DSRhoCPFit/DSRhoCPFit {params} --config={input.config} --plot-dir={output.plotdir} "
        "--output={output.result} &> {log}"

rule cpfit_data_systematics:
    input:
        config = "DSRhoCPFit/configs/{group}/numbered/config_data_{bkg}_{configno}.json",
    output:
        result = "DSRhoCPFit/results/{group}/{channel}_{type}_data_{bkg}/{configno}",
    log:
        "DSRhoCPFit/logs/{group}/{channel}_{type}_data_{bkg}/{configno}.log"
    params:
        "--MC=0",
        "--cpus=1",
        "--log",
        "--components=all",
        lambda wildcards:
            "--time-independent" if wildcards.type == "ti" else "",
        lambda wildcards:
            get_excluded_channels(wildcards)
    shell:
        "./DSRhoCPFit/DSRhoCPFit {params} --config={input.config} --output={output.result} &> {log}"

