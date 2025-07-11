#include "Input.hpp"
#include <random>

namespace boost_po = boost::program_options;

Input::Input(int argc, char *argv[])
{
	boost_po::options_description desc("Alowed options");
	desc.add_options()
		("help", "produce help message")
		("traj, j",     boost_po::value<std::string>(), "path to an input trajectory file")
		("toml, m",     boost_po::value<std::string>(), "path to a toml forcefield file")
		("output, o",   boost_po::value<std::string>(), "path to output")
		("cutoff, c",   boost_po::value<float>()->default_value(120.0),    "cutoff length")
		("bond, b",     boost_po::value<float>()->default_value(100.0),    "bond strength")
		("r0, r",       boost_po::value<float>()->default_value( 10.0),    "native bond length")
		("dihedral, d", boost_po::value<float>()->default_value(  1.0),    "dihedral strength")
		("theta0, t",   boost_po::value<float>()->default_value(  0.0),    "native torsion angle")
		("seed, s",     boost_po::value<int>()->default_value(-1),         "random seed")
		("maxcon, x",   boost_po::value<std::size_t>()->default_value(100),"max contact number");
	boost_po::variables_map vm;
	try
	{
		boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
	}
	catch(const boost_po::error_with_option_name& e)
	{
		std::cerr << e.what() << std::endl;
		std::exit(1);
	}
	boost_po::notify(vm);
	if (vm.count("help") || !vm.count("traj") || !vm.count("toml") || !vm.count("output"))
	{
		std::cerr << desc << std::endl;
		std::exit(1);
	}
	this->trajectory_name_ = vm["traj"].as<std::string>();
	this->base_ff_name_    = vm["toml"].as<std::string>();
	this->output_name_     = vm["output"].as<std::string>();
	this->cutoff_          = vm["cutoff"].as<float>();
	this->bond_k_          = vm["bond"].as<float>();
	this->r0_              = vm["r0"].as<float>();
	this->dihedral_k_      = vm["dihedral"].as<float>();
	this->theta0_          = vm["theta0"].as<float>();
	if (vm["seed"].as<int>() < 0)
	{
		std::random_device rng;
		this->seed_        = rng();
	}
	else this->seed_       = vm["seed"].as<int>();
	this->max_contact_     = vm["maxcon"].as<std::size_t>();
}
