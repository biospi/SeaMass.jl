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
      smbSpectrumIn.counts ./
        (smbSpectrumIn.locations[2:end] - smbSpectrumIn.locations[1:end-1]),
      0.0,
    ),
    line = :steppost,
    label = "input",
    title = filename * " SMB Comparison",
    xlabel = "m/z (Th)",
    ylabel = "ion count density",
    reuse = false,
  )
  plot!(
    smbSpectrumOutReconstruct.locations,
    vcat(
      smbSpectrumOutReconstruct.counts ./
        (smbSpectrumOutReconstruct.locations[2:end] - smbSpectrumOutReconstruct.locations[1:end-1]),
      0.0,
    ),
    line = :steppost,
    label = "seaMass-restore --reconstruct",
  )
  plot!(
    smbSpectrumOut.locations,
    smbSpectrumOut.counts,
    label = "seaMass-restore",
  )
  plot!(
    smbSpectrumOutDeconvolve.locations,
    smbSpectrumOutDeconvolve.counts,
    label = "seaMass-restore --deconvolve",
  )
  sticks!(
    smbSpectrumOutCentroid.locations,
    smbSpectrumOutCentroid.counts,
    label = "seaMass-restore --centroid",
    m = 4,
  )
  gui()

  # load mzMLb input spectrum
  mzmlbSpectrumIn = SeaMass.MzmlbSpectrum(
    "data/" * filename * ".mzMLb",
    inputSpectrumID
  )

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
    mzmlbSpectrumIn.mzs,
    mzmlbSpectrumIn.intensities,
    label = "input",
    title = "HYE124_TTOF6600_64var_lgillet_I150211_008__index_59994 mzMLb Comparison",
    xlabel = "m/z (Th)",
    ylabel = "intensity",
    reuse = false,
  )
  plot!(
    mzmlbSpectrumOut.mzs,
    mzmlbSpectrumOut.intensities,
    label = "seaMass-restore",
  )
  sticks!(
    mzmlbSpectrumOutCentroid.mzs,
    mzmlbSpectrumOutCentroid.intensities,
    label = "seaMass-restore --centroid",
    m = 4,
  )
  gui()

end

test("HYE124_TTOF6600_64var_lgillet_I150211_008__index_59994", "p-55-227433333333", 1, 1)
test("P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500", "p", (67 - 1) * 51 + 1, 67)

println("Press <Enter> to finish")
readline()
