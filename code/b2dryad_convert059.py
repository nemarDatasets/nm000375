"""IEEG059 (Dryad doi:10.5061/dryad.jdfn2z3k4; Combrisson et al. 2024 eLife) -> BIDS derivative dataset.

The release holds ONLY preprocessed high-gamma power (no raw iEEG): per subject a netCDF array
(n_trials, n_bipolar_contacts, 513 samples at 256 Hz, -0.5..1.5 s around outcome), float32, plus per-contact ROI labels
and per-trial behaviour (condition, prediction error, outcome). Representation: trials concatenated back-to-back in one
BrainVision IEEE_FLOAT_32 file per subject (values identical to the source float32, resolution 1, units n/a); one events
row per trial (onset = trial start in the concatenated file, plus the time-zero sample), carrying every behavioural
column verbatim. Nothing is filtered, resampled or rescaled.
Usage: python b2dryad_convert059.py <src_dir> <bids_root>
"""
import csv, glob, hashlib, io, json, os, re, shutil, sys, zipfile
import h5py, numpy as np, openpyxl

SRC, OUT = sys.argv[1], sys.argv[2]
Z = zipfile.ZipFile(glob.glob(os.path.join(SRC, "pblt_dataset*.zip"))[0])
os.makedirs(OUT, exist_ok=True)


def wtsv(p, h, rows):
    with open(p, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(h)
        for r in rows:
            w.writerow(["n/a" if (v is None or v == "") else v for v in r])


def wjson(p, o):
    with open(p, "w") as f:
        json.dump(o, f, indent=2, ensure_ascii=False)
        f.write("\n")


def xlsx(member):
    wb = openpyxl.load_workbook(io.BytesIO(Z.read(member)), read_only=True, data_only=True)
    rows = list(wb.worksheets[0].iter_rows(values_only=True))
    return [str(h) for h in rows[0]], rows[1:]


subs = sorted({int(re.search(r"power_subject-(\d+)\.nc", n).group(1)) for n in Z.namelist() if n.endswith(".nc")})
report = {"subjects": []}
parts = []
for s in subs:
    sub = f"sub-{s:02d}"
    d = os.path.join(OUT, sub, "ieeg")
    os.makedirs(d, exist_ok=True)
    f = h5py.File(io.BytesIO(Z.read(f"pblt_dataset/power/power_subject-{s}.nc")), "r")
    v = f["__xarray_dataarray_variable__"]
    x = v[()]
    attrs = {k: (a.decode() if isinstance(a, bytes) else (a.tolist() if hasattr(a, "tolist") else a)) for k, a in v.attrs.items() if not k.startswith("_") and k not in ("DIMENSION_LIST",)}
    ch = [c.decode() if isinstance(c, bytes) else str(c) for c in f["channels"][()]]
    times = f["times"][()]
    trials = f["trials"][()]
    assert x.dtype == np.float32 and x.ndim == 3
    ntr, nch, nt = x.shape
    sf = float(np.asarray(attrs["sfreq"]).ravel()[0])
    assert abs((times[1] - times[0]) * sf - 1) < 1e-9
    stem = f"{sub}_task-pblt_ieeg"
    cont = np.ascontiguousarray(x.transpose(0, 2, 1).reshape(ntr * nt, nch))  # (samples, channels) multiplexed
    cont.astype("<f4").tofile(os.path.join(d, stem + ".eeg"))
    vh = ["Brain Vision Data Exchange Header File Version 1.0", "; Written by b2dryad_convert059.py: authors' float32 high-gamma power, trials concatenated, resolution 1", "",
          "[Common Infos]", "Codepage=UTF-8", f"DataFile={stem}.eeg", f"MarkerFile={stem}.vmrk", "DataFormat=BINARY", "DataOrientation=MULTIPLEXED",
          f"NumberOfChannels={nch}", f"SamplingInterval={1e6 / sf!r}", "", "[Binary Infos]", "BinaryFormat=IEEE_FLOAT_32", "", "[Channel Infos]"]
    vh += [f"Ch{i}={c.replace(',', chr(92) + '1')},,1,n/a" for i, c in enumerate(ch, 1)]
    open(os.path.join(d, stem + ".vhdr"), "w", encoding="utf-8").write("\n".join(vh) + "\n")
    t0 = int(np.argmin(np.abs(times)))
    vm = ["Brain Vision Data Exchange Marker File, Version 1.0", "", "[Common Infos]", "Codepage=UTF-8", f"DataFile={stem}.eeg", "", "[Marker Infos]", "Mk1=New Segment,,1,1,0"]
    hb, rb = xlsx(f"pblt_dataset/behavior/subject-{s}.xlsx")
    assert len(rb) == ntr, (s, len(rb), ntr)
    ib = {h: i for i, h in enumerate(hb)}
    ev = []
    for k in range(ntr):
        r = rb[k]
        cond = r[ib["condition"]]
        ev.append(["%.8f" % (k * nt / sf), "%.8f" % (nt / sf), "trial", cond, k * nt, k + 1, int(trials[k]), r[ib["PE"]], r[ib["outcome_index"]], r[ib["outcome"]], "%.8f" % ((k * nt + t0) / sf), k * nt + t0])
        vm.append(f"Mk{len(vm) - 6}=Stimulus,{cond}_outcome,{k * nt + t0 + 1},1,0")
    open(os.path.join(d, stem + ".vmrk"), "w", encoding="utf-8").write("\n".join(vm) + "\n")
    wtsv(os.path.join(d, f"{sub}_task-pblt_events.tsv"),
         ["onset", "duration", "trial_type", "condition", "sample", "trial_index", "source_trial_label", "PE", "outcome_index", "outcome", "outcome_time_onset", "outcome_time_sample"], ev)
    ha, ra = xlsx(f"pblt_dataset/anatomy/subject-{s}.xlsx")
    ia = {h: i for i, h in enumerate(ha)}
    an = {r[ia["contacts"]]: (r[ia["roi"]], r[ia["hemisphere"]]) for r in ra}
    assert [r[ia["contacts"]] for r in ra] == ch, s
    wtsv(os.path.join(d, f"{sub}_task-pblt_channels.tsv"),
         ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "reference", "roi", "hemisphere", "status", "description"],
         [[c, "SEEG", "n/a", "n/a", "n/a", sf, "bipolar", an[c][0], an[c][1], "good", "high-gamma power of the bipolar derivation (authors' preprocessing)"] for c in ch])
    contacts = []
    for c in ch:
        for part in c.split("-"):
            if part not in contacts:
                contacts.append(part)
    wtsv(os.path.join(d, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "hemisphere"],
         [[c, "n/a", "n/a", "n/a", "n/a", re.sub(r"\d+$", "", c), next((an[b][1] for b in ch if c in b.split("-")), "n/a")] for c in contacts])
    wjson(os.path.join(d, f"{sub}_space-Other_coordsystem.json"), {
        "iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a",
        "iEEGCoordinateSystemDescription": "No electrode coordinates are included in the Dryad release; electrodes.tsv lists the contacts named in the bipolar channel labels (x, y, z = n/a). Regions of interest per bipolar derivation are in channels.tsv (column roi)."})
    wjson(os.path.join(d, stem + ".json"), {
        "TaskName": "pblt", "TaskDescription": "Probabilistic learning task with reward (gain) and punishment (loss) conditions; trials time-locked to the outcome (-0.5 to 1.5 s).",
        "SamplingFrequency": sf, "PowerLineFrequency": 50, "SoftwareFilters": "n/a", "HardwareFilters": "n/a",
        "iEEGReference": "bipolar (adjacent contacts on the same SEEG shaft, as named)",
        "RecordingType": "epoched", "EpochLength": nt / sf, "RecordingDuration": ntr * nt / sf,
        "SEEGChannelCount": nch,
        "DerivativeDescription": "Authors' high-gamma power (netCDF attributes: " + json.dumps(attrs) + "); trials concatenated back-to-back, values unchanged.",
    })
    report["subjects"].append({"sub": sub, "trials": ntr, "channels": nch, "samples_per_trial": nt, "sfreq": sf, "t0_index": t0,
                               "sha256_eeg": hashlib.sha256(open(os.path.join(d, stem + ".eeg"), "rb").read()).hexdigest(), "attrs": attrs,
                               "trial_labels": sorted(set(int(t) for t in trials))})
    parts.append([sub, f"subject-{s}", ntr, nch])
    print(sub, ntr, nch, nt, sf, flush=True)
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "source_id", "n_trials", "n_bipolar_contacts", "age", "sex", "handedness"],
     [p + ["n/a", "n/a", "n/a"] for p in parts])
wjson(os.path.join(OUT, "participants.json"), {"source_id": {"Description": "Subject label in the release file names"},
                                               "n_trials": {"Description": "Trials in the released high-gamma array"},
                                               "n_bipolar_contacts": {"Description": "Bipolar derivations in the released array"},
                                               "age": {"Description": "not given in the Dryad release", "Units": "year"}, "sex": {"Description": "not given in the Dryad release"},
                                               "handedness": {"Description": "not given in the Dryad release"}})
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
shutil.copy(__file__, os.path.join(OUT, "code", os.path.basename(__file__)))
json.dump(report, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
dst = os.path.join(OUT, "sourcedata", "dryad-jdfn2z3k4")
os.makedirs(dst, exist_ok=True)
for fn in os.listdir(SRC):
    p = os.path.join(SRC, fn)
    if os.path.isfile(p):
        shutil.copy(p, dst)
shutil.copytree(os.path.join(SRC, "_repository_metadata"), os.path.join(dst, "_repository_metadata"), dirs_exist_ok=True)
print("DONE", len(subs))
