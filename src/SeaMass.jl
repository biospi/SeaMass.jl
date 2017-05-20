module SeaMass


import HDF5
import LibExpat
using PyCall

const scipy_int = PyNULL()

function __init__()
  copy!(scipy_int, pyimport_conda("scipy.interpolate", "scipy"))
end


export MzmlSpectrum, SmiSpectrum, SmvSpectrum, SmoSpectrum, BinnedSmoSpectrum


type MzmlSpectrum
  mzs
  intensities
end

function MzmlSpectrum(filename::AbstractString, spectrumID::Number=1)
  HDF5.h5open(filename, "r") do file
    spectrumIndex = file["mzML_spectrumIndex"][spectrumID:spectrumID+1] + 1

    mzML = LibExpat.xp_parse(String(file["mzML"][spectrumIndex[1]:spectrumIndex[2]-1]))

    arrayLength = parse(Int, mzML["@defaultArrayLength"][1])
    if (arrayLength > 0)
      mzsDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@externalDataset"][1]
      mzsOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@offset"][1]) + 1
      intensitiesDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@externalDataset"][1]
      intensitiesOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@offset"][1]) + 1

      MzmlSpectrum(
        file[mzsDataset][mzsOffset:mzsOffset+arrayLength-1],
        file[intensitiesDataset][intensitiesOffset:intensitiesOffset+arrayLength-1]
      )
    else
      MzmlSpectrum(Float64[], Float32[])
    end
  end
end


type SmbSpectrum
  locations
  counts
  exposure
end

function SmbSpectrum(filename::AbstractString, spectrumID::Number=1)
  HDF5.h5open(filename, "r") do file
    if (HDF5.has(file,"countsIndex"))
      countsIndex = file["countsIndex"][spectrumID:spectrumID+1] + 1
    else
      countsIndex = [1 size(file["counts"])[1]+1]
    end

    if (countsIndex[2] - countsIndex[1] > 0)
      if (HDF5.has(file,"binLocations"))
        locations = file["binLocations"][countsIndex[1]+spectrumID-1:countsIndex[2]+spectrumID-1]
      else
        if(HDF5.has(file,"sampleLocations"))
          locations = file["sampleLocations"][countsIndex[1]:countsIndex[2]-1]
        else
          locations = file["centroidLocations"][countsIndex[1]:countsIndex[2]-1]
        end
      end

      if (HDF5.has(file,"exposures"))
        exposures = file["exposures"][spectrumID][1]
      else
        exposures = 1.0
      end

      SmbSpectrum(
        locations,
        file["counts"][countsIndex[1]:countsIndex[2]-1],
        exposures
      )         
    else
       SmbSpectrum(Float64[], Float32[], Float32[])  
    end
  end
end


type SmvSpectrum
    residualBinCounts
end

function SmvSpectrum(filename::AbstractString, spectrumID::Number=1)
  HDF5.h5open(filename, "r") do file
    if (HDF5.has(file,"binCountsIndex"))
      binCountsIndex = file["binCountsIndex"][spectrumID:spectrumID+1] + 1
    else
      binCountsIndex = [1 size(file["binCounts"])[1]+1]
    end

    SmvSpectrum(
      file["binCounts"][binCountsIndex[1]:binCountsIndex[2]-1],
    )
  end
end


type SmoSpectrum
  controlPoints
  offset
  scale
end

function SmoSpectrum(filename::AbstractString, spectrumID::Number=1)
  HDF5.h5open(filename, "r") do file
    dset = HDF5.d_open(file, "controlPoints")
    SmoSpectrum(
      transpose(file["controlPoints"][:,spectrumID])[1,:],
      HDF5.a_read(dset, "offset")[1],
      HDF5.a_read(dset, "scale")[1]
    )
    end
end


type BinnedSmoSpectrum
  binEdges
  binCounts
  binWidth
end

function BinnedSmoSpectrum(smoSpectrum::SmoSpectrum)
  binWidth = 1.0033548378 / (60 * 2^smoSpectrum.scale)
  binRange = (smoSpectrum.offset + [0, length(smoSpectrum.controlPoints) - 3]) * binWidth

  BinnedSmoSpectrum(
    linspace(binRange[1], binRange[2], length(smoSpectrum.controlPoints) - 2),
    filt([0.04f0, 0.46f0, 0.46f0, 0.04f0], 1, smoSpectrum.controlPoints)[4:end],
    binWidth
  )
end


function MzmlSpectrum(smoSpectrum::SmoSpectrum, nPoints::Number)
  coefs = vcat(smoSpectrum.controlPoints, 0, 0, 0, 0)
  knots = -3:length(coefs)-4
  tck = (knots, coefs, 3)
  x = linspace(0, length(smoSpectrum.controlPoints)-3, nPoints)

  binWidth = 1.0033548378 / (60 * 2^smoSpectrum.scale)
  mzRange = (smoSpectrum.offset + [0, length(smoSpectrum.controlPoints) - 3]) * binWidth

  MzmlSpectrum(
    linspace(mzRange[1], mzRange[2], nPoints),
    scipy_int[:splev](x, tck, ext=1) ./ binWidth
  )
end


end # module
