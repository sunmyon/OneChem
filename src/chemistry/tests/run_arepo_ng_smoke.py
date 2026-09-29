#!/usr/bin/env python3
"""NG integration smoke: hydro/EOS, snapshot, restart, and network validation.
Requires numpy, h5py, mpiexec and an NG CHEMISTRY+HELIUM+HLLD binary.
Usage: python3 run_arepo_ng_smoke.py /absolute/Arepo2 /absolute/fresh-output-dir
"""
import pathlib, shutil, subprocess, sys
import numpy as np
import h5py

executable = pathlib.Path(sys.argv[1]).resolve()
root = pathlib.Path(sys.argv[2]).resolve()
root.mkdir(parents=True, exist_ok=False)
chemistry = pathlib.Path(__file__).resolve().parents[1]
parameters = dict(line.split() for line in """InitCondFile /tmp/arepo_ng_chemistry_smoke/ics
OutputDir /tmp/arepo_ng_chemistry_smoke/output
SnapshotFileBase snap
OutputListFilename ol
ICFormat 3
SnapFormat 3
TimeLimitCPU 60
CpuTimeBetRestartFile 1800
MaxMemSize 1024
TopNodeFactor 4
ActivePartFracForNewDomainDecomp 0.01
TimeBegin 0.0
TimeMax 0.002
BoxSize 1.0
ComovingIntegrationOn 0
Omega0 0
OmegaLambda 0
OmegaBaryon 0
OmegaCurvature 0
Hubble 0
HubbleParam 1.0
OutputListOn 0
TimeBetSnapshot 0.001
TimeOfFirstSnapshot 0.0
TimeBetStatistics 0.001
NumFilesPerSnapshot 1
MaxFilesWithConcurrentIO 1
ErrTolIntAccuracy 0.025
CourantFac 0.4
MaxSizeTimestep 0.001
MinSizeTimestep 0
InitGasTemp 1000
MinGasTemp 0
MinEgySpec 1e-6
TypeOfOpeningCriterion 1
ErrTolTheta 0.7
ErrTolThetaMax 0.9
ErrTolForceAcc 0.0025
NumberOfGlobalTimestepsToStopAfter 0
DesNumNgb 64
MaxNumNgbDeviation 2
UnitLength_in_cm 1e12
UnitMass_in_g 1e16
UnitVelocity_in_cm_per_s 1e5
GravityConstantInternal 1e-20
SofteningComovingClass0 7.81250e-05
SofteningComovingClass1 0.0
SofteningMaxPhysClass0 7.81250e-05
SofteningMaxPhysClass1 0.0
SofteningClassOfPartType0 0
SofteningClassOfPartType1 1
SofteningClassOfPartType2 1
SofteningClassOfPartType3 1
SofteningClassOfPartType4 1
SofteningClassOfPartType5 1
CellMaxAngleFactor 2
CellShapingSpeed 0.25
MinimumComovingHydroSoftening 0.001
AdaptiveHydroSofteningSpacing 1.2
GasSoftFactor 2.5
ChemistryTableRoot /Users/j/Desktop/codex/1ZONE/src/chemistry
ChemistryInitialMetallicity 0
ChemistryInitialIonizedFraction 1e-4
ChemistryInitialMolecularFraction 0.1
ChemistryCMBTemperature 2.728
ChemistryBalanceWork 1
""".strip().splitlines())
parameters.update(InitCondFile=str(root / "ics"), OutputDir=str(root / "output"), ChemistryTableRoot=str(chemistry))

def parameter_file(name, **changes):
    p = parameters | changes
    path = root / (name + ".param")
    path.write_text("".join(f"{key} {value}\n" for key, value in p.items()))
    return path

def run(name, param, *args, expected_failure=None):
    with (root / (name + ".log")).open("w") as log:
        result = subprocess.run(["mpiexec", "-n", "1", str(executable), str(param), *args], cwd=root, stdout=log, stderr=subprocess.STDOUT, timeout=120)
    text = (root / (name + ".log")).read_text(errors="replace")
    if expected_failure:
        assert result.returncode != 0 and expected_failure in text, name
    else:
        assert result.returncode == 0 and "Final time=" in text, name

n = 8
xyz = np.array(np.meshgrid(*([(np.arange(n)+0.5)/n]*3), indexing="ij")).reshape(3,-1).T
count = len(xyz)
with h5py.File(root / "ics.hdf5", "w") as f:
    h = f.create_group("Header")
    counts = np.array([count,0,0,0,0,0], dtype=np.uint32)
    for key, value in dict(NumPart_ThisFile=counts, NumPart_Total=counts,
        NumPart_Total_HighWord=np.zeros(6,dtype=np.uint32), MassTable=np.zeros(6),
        Time=0., Redshift=0., BoxSize=1., NumFilesPerSnapshot=1, Omega0=0.,
        OmegaLambda=0., HubbleParam=1., Flag_DoublePrecision=1, Flag_Entropy_ICs=0).items():
        h.attrs[key] = value
    gas = f.create_group("PartType0")
    gas["Coordinates"] = xyz
    gas["Velocities"] = np.zeros((count,3))
    gas["ParticleIDs"] = np.arange(1,count+1,dtype=np.uint64)
    gas["Masses"] = np.full(count,1./count)
    gas["InternalEnergy"] = np.full(count,10.)

run("initial", parameter_file("initial"))
run("restart", parameter_file("restart", TimeMax="0.004"), "1")
run("continuous", parameter_file("continuous", TimeMax="0.004", OutputDir=str(root / "continuous")))
run("snapshot", parameter_file("snapshot", TimeMax="0.004", OutputDir=str(root / "snapshot"), InitCondFile=str(root / "output/snap_002")), "2")

fields = ["ChemicalAbundances", "ChemistryGamma", "ChemistryMetallicity", "InternalEnergy", "Pressure", "Masses"]
with h5py.File(root / "continuous/snap_004.hdf5") as baseline:
    gas = baseline["PartType0"]
    assert np.isfinite(gas["ChemicalAbundances"][:]).all()
    assert (gas["ChemicalAbundances"][:] >= 0).all()
    assert (gas["ChemistryGamma"][:] < 5./3.).all()
    assert np.allclose(gas["Pressure"][:], (gas["ChemistryGamma"][:]-1)*gas["Density"][:]*gas["InternalEnergy"][:], rtol=1e-12)
    assert np.allclose(gas["Masses"][:].sum(), 1., rtol=1e-12)
    for folder in ["output", "snapshot"]:
        with h5py.File(root / folder / "snap_004.hdf5") as resumed:
            for field in fields:
                np.testing.assert_allclose(gas[field][:], resumed["PartType0"][field][:], rtol=1e-12, atol=1e-50)
    with h5py.File(root / "output/snap_000.hdf5") as initial:
        assert np.any(gas["InternalEnergy"][:] != initial["PartType0/InternalEnergy"][:])

shutil.copyfile(root / "output/snap_002.hdf5", root / "incompatible.hdf5")
with h5py.File(root / "incompatible.hdf5", "r+") as f:
    d = f["PartType0/ChemistryNetworkSignature"]
    d[...] = d[:] ^ np.uint64(1)
run("incompatible", parameter_file("incompatible", TimeMax="0.004", OutputDir=str(root / "rejected"), InitCondFile=str(root / "incompatible")), "2", expected_failure="network signature")
# Exercise nonzero hydro fluxes with spatially varying molecular fraction and Z.
shutil.copyfile(root / "output/snap_000.hdf5", root / "inhomogeneous.hdf5")
with h5py.File(root / "inhomogeneous.hdf5", "r+") as f:
    g = f["PartType0"]
    left = g["Coordinates"][:,0] < 0.5
    a = g["ChemicalAbundances"][:]
    assert a.shape[1] == 8, "This fixture requires precisely the H/He network"
    a[left,6] = 0.01  # H2, in the H/He species ordering
    a[:,1] = 0.76 - a[:,2] - 2*a[:,6] - a[:,7]
    g["ChemicalAbundances"][...] = a
    g["ChemistryMetallicity"][...] = left.astype(float)
    g["InternalEnergy"][...] = np.where(left, 15., 10.)
run("transport", parameter_file("transport", OutputDir=str(root / "transport"), InitCondFile=str(root / "inhomogeneous")))
with h5py.File(root / "transport/snap_000.hdf5") as initial, h5py.File(root / "transport/snap_002.hdf5") as final:
    before, after = initial["PartType0"], final["PartType0"]
    mass0, mass1 = before["Masses"][:], after["Masses"][:]
    assert np.max(np.abs(np.sort(mass1)-np.sort(mass0))) > 1.e-9, "No hydro mass transport occurred"
    assert (after["ChemicalAbundances"][:] >= 0).all()
    np.testing.assert_allclose(mass1.sum(), mass0.sum(), rtol=1.e-12)
    np.testing.assert_allclose(np.dot(mass1, after["ChemistryMetallicity"][:]), np.dot(mass0, before["ChemistryMetallicity"][:]), rtol=1.e-12)
    def hydrogen(g):
        a = g["ChemicalAbundances"][:]
        return np.dot(g["Masses"][:], a[:,1]+a[:,2]+2*a[:,6]+a[:,7])
    np.testing.assert_allclose(hydrogen(after), hydrogen(before), rtol=1.e-10)

# Missing thermal energy must use the composition-dependent temperature floor.
with h5py.File(root / "ics.hdf5", "r+") as f:
    f["PartType0/InternalEnergy"][...] = 0.
run("temperature_floor", parameter_file("temperature_floor", OutputDir=str(root / "temperature_floor"), InitGasTemp="50", MinGasTemp="500"))
with h5py.File(root / "temperature_floor/snap_000.hdf5") as f:
    g = f["PartType0"]
    mu = 1. / g["ChemicalAbundances"][:].sum(axis=1)
    temperature = g["InternalEnergy"][:] * 1.e10 * (g["ChemistryGamma"][:]-1) * 1.6726e-24 * mu / 1.3806e-16
    np.testing.assert_allclose(temperature, 500., rtol=1.e-12)
print("NG chemistry: EOS, evolving abundances, snapshot and binary restart equivalence, incompatible-network rejection, conservative transport, temperature floor passed")
