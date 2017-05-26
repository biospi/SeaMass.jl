module SeaMass


import HDF5
import LibExpat


export MzmlbSpectrum, SmbSpectrum


type MzmlbSpectrum
  mzs
  intensities
end

function MzmlbSpectrum(filename::AbstractString, spectrumID::Number=1)
  HDF5.h5open(filename, "r") do file
    spectrumIndex = file["mzML_spectrumIndex"][spectrumID:spectrumID+1] + 1

    mzML = LibExpat.xp_parse(String(file["mzML"][spectrumIndex[1]:spectrumIndex[2]-1]))

    arrayLength = parse(Int, mzML["@defaultArrayLength"][1])
    if (arrayLength > 0)
      mzsDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@externalDataset"][1]
      mzsOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@offset"][1]) + 1
      intensitiesDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@externalDataset"][1]
      intensitiesOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@offset"][1]) + 1

      MzmlbSpectrum(
        file[mzsDataset][mzsOffset:mzsOffset+arrayLength-1],
        file[intensitiesDataset][intensitiesOffset:intensitiesOffset+arrayLength-1]
      )
    else
      MzmlbSpectrum(Float64[], Float32[])
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


end # module
