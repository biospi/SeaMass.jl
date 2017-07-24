cd(ENV["PWD"])

import SeaMass
using Base.Test
using Plots
pyplot()

function test(filename, id, inputSpectrumID, outputSpectrumID)

  # load SMB input spectrum
  smbSpectrumIn = SeaMass.SmbSpectrum(
    "data/out/" * filename * "/1.seamass/" * filename * "." * id * ".smb",
    outputSpectrumID
  )

  # load SMB output spectrum
  smbSpectrumOut = SeaMass.SmbSpectrum(
    "data/out/" * filename * "/2.seamass-restore/" * filename * "." * id * ".smb",
    outputSpectrumID
  )

  # load SMB deconvolved output spectrum
  smbSpectrumOutDeconvolve = SeaMass.SmbSpectrum(
    "data/out/" * filename * "/3.seamass-restore_--deconvolve/" * filename * "." * id * ".smb",
    outputSpectrumID
  )

  # load SMB reconstructed output spectrum
  smbSpectrumOutReconstruct = SeaMass.SmbSpectrum(
    "data/out/" * filename * "/4.seamass-restore_--reconstruct/" * filename * "." * id * ".smb",
    outputSpectrumID
  )

  # load SMB centroided output spectrum
  smbSpectrumOutCentroid = SeaMass.SmbSpectrum(
    "data/out/" * filename * "/5.seamass-restore_--centroid/" * filename * "." * id * ".smb",
    outputSpectrumID
  )

  # Compare smb
  plot(
    smbSpectrumIn.locations,
    vcat(
      smbSpectrumIn.counts ./ (smbSpectrumIn.locations[2:end] - smbSpectrumIn.locations[1:end-1]),
      0.0,
    ),
    line = :steppost,
    label = "SMB input",
    title = "SMB comparison - " * filename,
    xlabel = "m/z (Th)",
    ylabel = "ion count density",
    reuse = false,
  )
  plot!(
    smbSpectrumOutReconstruct.locations,
    vcat(
      smbSpectrumOutReconstruct.counts ./  (smbSpectrumOutReconstruct.locations[2:end] - smbSpectrumOutReconstruct.locations[1:end-1]),
      0.0,
    ),
    line = :steppost,
    label = "SMB output (seaMass-restore --reconstruct)",
  )
  plot!(
    smbSpectrumOut.locations,
    smbSpectrumOut.counts[2:end],
    label = "SMB output (seaMass-restore)",
  )
  plot!(
    smbSpectrumOutDeconvolve.locations,
    smbSpectrumOutDeconvolve.counts,
    label = "SMB output (seaMass-restore --deconvolve)",
  )
  sticks!(
    smbSpectrumOutCentroid.locations,
    smbSpectrumOutCentroid.counts,
    label = "SMB output (seaMass-restore --centroid)",
    m = 4,
  )
  gui()

  # load mzMLb output spectrum
  mzmlbSpectrumOut = SeaMass.MzmlbSpectrum(
    "data/out/" * filename * "/2.seamass-restore/" * filename * ".mzMLb",
    outputSpectrumID
  )

  # load mzMLb centroided output spectrum
  mzmlbSpectrumOutCentroid = SeaMass.MzmlbSpectrum(
    "data/out/" * filename * "/5.seamass-restore_--centroid/" * filename * ".mzMLb",
    outputSpectrumID
  )

  # Compare mzMLb
  plot(
    smbSpectrumIn.locations,
    vcat(
      smbSpectrumIn.counts ./ (smbSpectrumIn.locations[2:end] - smbSpectrumIn.locations[1:end-1]),
      0.0,
    ),
    line = :steppost,
    label = "SMB input",
    title = "mzMLb comparison - " * filename,
    xlabel = "m/z (Th)",
    ylabel = "intensity density",
    reuse = false,
  )
  plot!(
    mzmlbSpectrumOut.mzs,
    mzmlbSpectrumOut.intensities,
    label = "mzMLb output (seaMass-restore)",
  )
  sticks!(
    mzmlbSpectrumOutCentroid.mzs,
    mzmlbSpectrumOutCentroid.intensities,
    label = "mzMLb output (seaMass-restore --centroid)",
    m = 4,
  )
  gui()

end

filename = "HYE124_TTOF6600_64var_lgillet_I150211_008__index_59994"
id = "p-55-227433333333"
inputSpectrumID = 1
outputSpectrumID = 1
test(filename, id, inputSpectrumID, outputSpectrumID)

filename = "P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500"
id = "p"
inputSpectrumID = (67 - 1) * 51 + 1
outputSpectrumID = 67
test(filename, id, inputSpectrumID, outputSpectrumID)

println("Press <Enter> to finish")
readline()
