#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "NTK_aligner_api.hpp"
#include "NT_aligner_api.hpp"

namespace py = pybind11;
using dynamont::PoreType;

namespace
{

PoreType poreTypeFromString(const std::string& pore)
{
    if (pore == "rna002")
        return PoreType::RNA002;
    if (pore == "rna004")
        return PoreType::RNA004;
    if (pore == "dna_r9")
        return PoreType::DNA_R9;
    if (pore == "dna_r10_260bps")
        return PoreType::DNA_R10_260;
    if (pore == "dna_r10_400bps")
        return PoreType::DNA_R10_400;

    throw std::invalid_argument("Unknown pore type: " + pore);
}

std::unique_ptr<dynamont::Aligner> makeAligner(
    const std::string& modelFile,
    PoreType pore,
    const std::string& mode,
    int threads,
    std::size_t band)
{
    if (mode == "basic" || mode == "nt")
    {
        return std::make_unique<dynamont::NTAligner>(modelFile, pore, threads, band);
    }
    if (mode == "resquiggle" || mode == "ntk")
    {
        return std::make_unique<dynamont::NTKAligner>(modelFile, pore, threads, band);
    }

    throw std::invalid_argument("Unknown aligner mode: " + mode);
}

py::dict resultToPython(const dynamont::Result& result)
{
    py::dict output;
    output["Z"] = result.Z;

    py::array_t<std::size_t> sequencePositions(result.segments.size());
    py::array_t<std::size_t> signalPositions(result.segments.size());
    py::array_t<double> probabilities(result.segments.size());
    py::list states;
    py::list polishes;

    auto sequenceView = sequencePositions.mutable_unchecked<1>();
    auto signalView = signalPositions.mutable_unchecked<1>();
    auto probabilityView = probabilities.mutable_unchecked<1>();

    for (std::size_t i = 0; i < result.segments.size(); ++i)
    {
        const auto& segment = result.segments[i];
        sequenceView(i) = segment.sequencePosition;
        signalView(i) = segment.signalPosition;
        probabilityView(i) = segment.probability;
        states.append(py::str(std::string(1, segment.state)));
        polishes.append(py::str(segment.polish));
    }

    output["sequence_positions"] = std::move(sequencePositions);
    output["signal_positions"] = std::move(signalPositions);
    output["probabilities"] = std::move(probabilities);
    output["states"] = std::move(states);
    output["polishes"] = std::move(polishes);
    return output;
}

py::dict trainingResultToPython(const dynamont::TrainingResult& result)
{
    py::dict transitions;
    transitions["m1"] = result.transitions.m1;
    transitions["e1"] = result.transitions.e1;
    transitions["e2"] = result.transitions.e2;

    py::list emissionModel;
    for (const auto& entry : result.emissionModel)
    {
        py::dict item;
        item["mean"] = entry.mean;
        item["stdev"] = entry.stdev;
        emissionModel.append(std::move(item));
    }

    py::dict output;
    output["Z"] = result.Z;
    output["transition_params"] = std::move(transitions);
    output["emission_model"] = std::move(emissionModel);
    return output;
}

class PyAligner
{
public:
    PyAligner(
        const std::string& modelFile,
        PoreType pore,
        const std::string& mode,
        int threads,
        std::size_t band):
        aligner_(makeAligner(modelFile, pore, mode, threads, band))
    {
    }

    PyAligner(
        const std::string& modelFile,
        const std::string& pore,
        const std::string& mode,
        int threads,
        std::size_t band):
        PyAligner(modelFile, poreTypeFromString(pore), mode, threads, band)
    {
    }

    py::dict align(
        const py::array_t<double, py::array::c_style | py::array::forcecast>& signal,
        const std::string& sequence,
        bool calcProbabilities)
    {
        if (signal.ndim() != 1)
            throw std::invalid_argument("Signal must be a one-dimensional array");

        const auto signalLength = static_cast<std::size_t>(signal.shape(0));
        dynamont::Result result;
        {
            py::gil_scoped_release release;
            result = aligner_->align(signal.data(), signalLength, sequence, calcProbabilities);
        }
        return resultToPython(result);
    }

    py::dict train(
        const py::array_t<double, py::array::c_style | py::array::forcecast>& signal,
        const std::string& sequence)
    {
        if (signal.ndim() != 1)
            throw std::invalid_argument("Signal must be a one-dimensional array");

        const auto signalLength = static_cast<std::size_t>(signal.shape(0));
        dynamont::TrainingResult result;
        {
            py::gil_scoped_release release;
            result = aligner_->train(signal.data(), signalLength, sequence);
        }
        return trainingResultToPython(result);
    }

private:
    std::unique_ptr<dynamont::Aligner> aligner_;
};

py::dict align(
    PyAligner& aligner,
    const py::array_t<double, py::array::c_style | py::array::forcecast>& signal,
    const std::string& sequence,
    bool calcProbabilities)
{
    return aligner.align(signal, sequence, calcProbabilities);
}

} // namespace

PYBIND11_MODULE(_dynamont, module)
{
    module.doc() = "Native dynamont alignment API.";

    py::enum_<PoreType>(module, "PoreType")
        .value("RNA002", PoreType::RNA002)
        .value("RNA004", PoreType::RNA004)
        .value("DNA_R9", PoreType::DNA_R9)
        .value("DNA_R10_260", PoreType::DNA_R10_260)
        .value("DNA_R10_400", PoreType::DNA_R10_400);

    py::class_<PyAligner>(module, "Aligner")
        .def(
            py::init<const std::string&, PoreType, const std::string&, int, std::size_t>(),
            py::arg("model_file"),
            py::arg("pore"),
            py::arg("mode") = "basic",
            py::arg("threads") = 1,
            py::arg("band") = 400)
        .def(
            py::init<const std::string&, const std::string&, const std::string&, int, std::size_t>(),
            py::arg("model_file"),
            py::arg("pore"),
            py::arg("mode") = "basic",
            py::arg("threads") = 1,
            py::arg("band") = 400)
        .def(
            "align",
            &align,
            py::arg("signal"),
            py::arg("sequence"),
            py::arg("calc_probabilities") = false)
        .def(
            "train",
            &PyAligner::train,
            py::arg("signal"),
            py::arg("sequence"));

    module.def("pore_type", &poreTypeFromString, py::arg("pore"));
}
