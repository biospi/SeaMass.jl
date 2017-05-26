import SeaMass
using Base.Test
using Plots
pyplot()

cd(ENV["PWD"])

# spectrum ID of spectrum we will be plotting
spectrumID = 67

# load mzMLb input spectrum
mzmlbSpectrumIn = SeaMass.MzmlbSpectrum(
  "data/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.mzMLb",
  (spectrumID - 1) * 51 + 1
)

# load SMB input spectrum
smbSpectrumIn = SeaMass.SmbSpectrum(
  "data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500" *
    "/1.mzmlb2smb/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.p.smb",
  spectrumID
)

# load SMB output spectrum
smbSpectrumOut = SeaMass.SmbSpectrum(
  "data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500" *
    "/3.seamass-restore/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.p.smb",
  spectrumID
)

# load SMB centroided output spectrum
smbSpectrumOutCentroid = SeaMass.SmbSpectrum(
  "data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500" *
    "/5.seamass-restore_--centroid/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.p.smb",
  spectrumID
)

# load mzMLb output spectrum
mzmlbSpectrumOut = SeaMass.MzmlbSpectrum(
  "data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500" *
    "/4.smb2mzmlb/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.mzMLb",
  spectrumID
)

# load mzMLb output spectrum
mzmlbSpectrumOutCentroid = SeaMass.MzmlbSpectrum(
  "data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500" *
    "/6.smb2mzmlb/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.mzMLb",
  spectrumID
)

# Compare mzMLb input to 'mzmlb2smb' smb input (binned)
plot(
  mzmlbSpectrumIn.mzs,
  mzmlbSpectrumIn.intensities,
  m = 1,
  label = "mzMLb input",
  title = "Compare mzMLb input (sampled) to 'mzmlb2smb' smb input (binned)",
  xlabel = "m/z (Th)",
  ylabel = "intensity",
  xlims = (602, 605),
  reuse = false,
)
plot!(
  smbSpectrumIn.locations,
  vcat(smbSpectrumIn.counts / smbSpectrumIn.exposure, 0.0),
  line = :steppost,
  label = "SMB input",
)
gui()

# Compare smb input to 'seamass' ⇨ 'seamass-restore' smb output
plot(
  smbSpectrumIn.locations,
  vcat(smbSpectrumIn.counts, 0.0),
  line = :steppost,
  label = "SMB input",
  m = 1,
  title = "Compare smb input to 'seamass' ⇨ 'seamass-restore' smb output)",
  xlabel = "m/z (Th)",
  ylabel = "ion count",
  xlims = (602, 605),
  reuse = false,
)
plot!(
  smbSpectrumOut.locations,
  vcat(smbSpectrumOut.counts, 0.0),
  line = :steppost,
  label = "seaMass-restore SMB output",
)
gui()

# Compare smb output to 'seamass-restore' mzMLb output
plot(
  smbSpectrumOut.locations,
  vcat(
    smbSpectrumOut.counts ./
      (smbSpectrumOut.locations[2:end] - smbSpectrumOut.locations[1:end-1]),
    0.0,
  ),
  line = :steppost,
  label = "seaMass-restore SMB output",
  m = 1,
  title = "Compare smb output to 'seamass-restore' mzMLb output",
  xlabel = "m/z (Th)",
  ylabel = "ion count density",
  xlims = (602, 605),
  reuse = false,
)
plot!(
  mzmlbSpectrumOut.mzs,
  mzmlbSpectrumOut.intensities * smbSpectrumOut.exposure,
  m = 1,
  label = "mzMLb output",
)
gui()

# Compare smb input to 'seamass' ⇨ 'seamass-restore --centroid' smb output
plot(
  smbSpectrumIn.locations,
  vcat(smbSpectrumIn.counts, 0.0),
  line = :steppost,
  label = "SMB input",
  m = 1,
  title = "Compare smb input to 'seamass' ⇨ 'seamass-restore --centroid' smb output",
  xlabel = "m/z (Th)",
  ylabel = "ion count",
  xlims = (602, 605),
  reuse = false,
)
sticks!(
  smbSpectrumOutCentroid.locations,
  smbSpectrumOutCentroid.counts,
  label = "seaMass-restore SMB output",
  m = 4,
)
gui()

# Compare mzML input to 'seamass-restore --centroid' mzMLb output
plot(
  mzmlbSpectrumIn.mzs,
  mzmlbSpectrumIn.intensities * smbSpectrumOut.exposure,
  m = 1,
  label = "mzMLb input",
  m = 1,
  title = "Compare mzML input to 'seamass-restore --centroid' mzMLb output",
  xlabel = "m/z (Th)",
  ylabel = "ion count",
  xlims = (602, 605),
  reuse = false,
)
sticks!(
  smbSpectrumOutCentroid.locations,
  smbSpectrumOutCentroid.counts,
  label = "ion count density",
  m = 4,
)
gui()

println("Press <Enter> to finish")
readline()
