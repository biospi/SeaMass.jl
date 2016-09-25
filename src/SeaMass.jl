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

function MzmlSpectrum(filename::AbstractString, spectrumID::Number)
    HDF5.h5open(filename, "r") do file
        spectrumIndex = file["mzML_spectrumIndex"][spectrumID:spectrumID+1] + 1

        mzML = LibExpat.xp_parse(String(file["mzML"][spectrumIndex[1]:spectrumIndex[2]-1]))

        arrayLength = parse(Int, mzML["@defaultArrayLength"][1])
        mzsDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@externalDataset"][1]
        mzsOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@offset"][1]) + 1
        intensitiesDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@externalDataset"][1]
        intensitiesOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@offset"][1]) + 1

        MzmlSpectrum(
            file[mzsDataset][mzsOffset:mzsOffset+arrayLength-1],
            file[intensitiesDataset][intensitiesOffset:intensitiesOffset+arrayLength-1]
        )
    end
end


type SmiSpectrum
    binEdges
    binCounts
    exposure
end

function SmiSpectrum(filename::AbstractString, spectrumID::Number)
    HDF5.h5open(filename, "r") do file
        if (HDF5.has(file,"spectrumIndex"))
            spectrumIndex = file["spectrumIndex"][spectrumID:spectrumID+1] + 1
        else
            spectrumIndex = [1 size(file["binCounts"])[1]+1]
        end
        
        SmiSpectrum(
            file["binEdges"][spectrumIndex[1]+spectrumID-1:spectrumIndex[2]+spectrumID-1],
            file["binCounts"][spectrumIndex[1]:spectrumIndex[2]-1],
            file["exposures"][spectrumID][1]
        )
    end
end


type SmvSpectrum
    residualBinCounts
end

function SmvSpectrum(filename::AbstractString, spectrumID::Number)
    HDF5.h5open(filename, "r") do file
        if (HDF5.has(file,"spectrumIndex"))
            spectrumIndex = file["spectrumIndex"][spectrumID:spectrumID+1] + 1
        else
            spectrumIndex = [1 size(file["binCounts"])[1]+1]
        end
        
        SmvSpectrum(
            file["binCounts"][spectrumIndex[1]:spectrumIndex[2]-1],
        )
    end
end


type SmoSpectrum
  controlPoints
  offset
  scale
end

function SmoSpectrum(filename::AbstractString, spectrumID::Number)
    HDF5.h5open(filename, "r") do file
        SmoSpectrum(
            read(file, "controlPoints"),
            HDF5.h5readattr(filename, "controlPoints")["offset"][1],
            HDF5.h5readattr(filename, "controlPoints")["scale"][1]
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
