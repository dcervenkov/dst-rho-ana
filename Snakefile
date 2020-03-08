import glob

STREAMS = range(0, 6)
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
            else:
                return (MC_WO_SIGNAL["Kpi"] +
                        MC_WO_SIGNAL["Kpipi0"] +
                        MC_WO_SIGNAL["K3pi"])
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


rule all:
    input:
        yield_jobs = (
            expand("DSRhoYield/plots/{channel}_stream{stream}",
                   channel=CHANNELS, stream=STREAMS),
            expand("DSRhoYield/plots/{channel}", channel=CHANNELS)
            ),

        background_jobs = (
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
                   channel=CHANNELS, component=["all"], type=["lifetime", "mixing"],
                   stream=range(6)),
            expand("DSRhoLifetime/plots/{channel}_{component}_{type}_{stream}",
                   channel=CHANNELS, component=["all"], type=["lifetime", "mixing"],
                   stream=range(6)),

            expand("DSRhoLifetime/results/{channel}_{type}_data",
                   channel=CHANNELS_AND_TOGETHER, component=["all"], type=["lifetime", "mixing"])
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

rule yield_mc:
    input:
        "data/{channel}/realistic_mc",
        "data/{channel}/realistic_mc/stream{stream}"
    output:
        directory("DSRhoYield/plots/{channel}_stream{stream}")
    log:
        "DSRhoYield/log/{channel}_stream{stream}"
    shell:
        "./DSRhoYield/DSRhoYield --MC {input} {output} &> {log}"

rule yield_data:
    input:
        "data/{channel}/realistic_mc",
        "data/{channel}"
    output:
        directory("DSRhoYield/plots/{channel}")
    log:
        "DSRhoYield/log/{channel}"
    shell:
        "./DSRhoYield/DSRhoYield {input} {output} &> {log}"

rule yield_summary:
    input:
        expand("DSRhoYield/plots/{channel}_stream{stream}/fit_results.root",
               channel=CHANNELS, stream=STREAMS)
    output:
        expand("DSRhoYield/results/{channel}_{type}_fractions.json", channel=CHANNELS, type=["mc", "data"]),
        "DSRhoYield/results/avg_data_fractions.json"
    log:
        "yield.log"
    shell:
        "cd DSRhoYield && ./tools/print_all_yields.py &> {log}"

rule background:
    input:
        get_data_for_background
    output:
        result = "DSRhoBackground/results/{channel}_{mc}_{bkg_type}.json",
        plotdir = directory("DSRhoBackground/plots/{channel}_{mc}_{bkg_type}")
    log:
        "DSRhoBackground/logs/{channel}_{mc}_{bkg_type}"
    params:
        lambda wildcards:
            "--notail --nodelta --physics" if wildcards.bkg_type == "scf" else
            "--notail --physics"
    wildcard_constraints:
        channel = "Kpi|Kpipi0|K3pi|together"
    shell:
        "./DSRhoBackground/DSRhoBackground"
        " {params} --plot-dir={output.plotdir} {output.result} {input} &> {log}"

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
        " --plot-dir={output.plotdir} {output.result} {input} &> {log}"

rule lifetime_configs:
    input:
        expand(rules.background.output, channel=CHANNELS, mc="mc", bkg_type=BKG_TYPES),
        expand(rules.background.output, channel=CHANNELS, mc="data", bkg_type="sidebands"),
        rules.yield_summary.output,
        template = "DSRhoLifetime/configs/templates/{config}.template.json"
    output:
        "DSRhoLifetime/configs/{config}.json"
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

