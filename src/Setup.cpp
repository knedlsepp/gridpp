#include "Setup.h"
#include "File/File.h"
#include "Calibrator/Calibrator.h"
#include "Downscaler/Downscaler.h"

Setup::Setup(const std::vector<std::string>& argv) {

   std::string inputFilename = "";
   std::string outputFilename = "";

   // Process input/output filenames and options
   int index = 0;
   while(index < argv.size()) {
      std::string arg = argv[index];
      if(inputFilename == "") {
         inputFilename = arg;
      }
      else if(outputFilename == "") {
         if(Util::hasChar(arg, '=')) {
            inputOptions.addOptions(arg);
         }
         else {
            outputFilename = arg;
         }
      }
      else {
         if(Util::hasChar(arg, '=')) {
            outputOptions.addOptions(arg);
         }
         else {
            break;
         }
      }
      index++;
   }
   std::vector<std::string> inputFilenames = Util::glob(inputFilename);
   std::vector<std::string> outputFilenames = Util::glob(outputFilename);
   if(inputFilenames.size() != outputFilenames.size()) {
      std::stringstream ss;
      ss << "Unequal number of input (" << inputFilenames.size() << ") and output ("
         << outputFilenames.size() << ")";
      Util::error(ss.str());
   }
   if(inputFilenames.size() == 0) {
      Util::error("No valid input files");
   }
   if(outputFilenames.size() == 0) {
      Util::error("No valid output files");
   }

   for(int f = 0; f < outputFilenames.size(); f++) {
      if(!hasFile(outputFilenames[f])) {
         std::unique_ptr<File> outputFile = File::getScheme(outputFilenames[f], outputOptions, false);
         if(!outputFile) {
            Util::error("File '" + outputFilenames[f] + " is invalid");
         }
         mFileMap[outputFilenames[f]] = std::move(outputFile);
      }
      outputFiles.push_back(mFileMap[outputFilenames[f]].get());

      if(!hasFile(inputFilenames[f])) {
         std::unique_ptr<File> inputFile = File::getScheme(inputFilenames[f], inputOptions, true);
         if(inputFile == NULL) {
            Util::error("File '" + inputFilenames[f] + " is invalid");
         }
         mFileMap[inputFilenames[f]] = std::move(inputFile);
      }
      inputFiles.push_back(mFileMap[inputFilenames[f]].get());
   }

   // Implement a finite state machine
   enum State {START = 0, VAR = 1, VAROPT = 2, NEWVAR = 3, DOWN = 10, DOWNOPT = 15, CAL = 20, NEWCAL = 22, CALOPT = 25, PARDOWN = 30, PAROPTDOWN = 35, PARCAL = 40, PAROPTCAL = 45, END = 90, ERROR = 100};
   State state = START;
   State prevState = START;

   Variable::Type variable = Variable::None;
   Options vOptions;
   Options dOptions;
   Options cOptions;
   Options pOptions;
   std::string errorMessage = "";

   std::string downscaler = defaultDownscaler();
   std::string calibrator = "";
   std::string parameterFile = "";
   std::vector<std::unique_ptr<Calibrator>> calibrators;
   std::vector<std::unique_ptr<ParameterFile>> parameterFileCalibrators;
   std::unique_ptr<ParameterFile> parameterFileDownscaler;
   while(true) {
      // std::cout << state << std::endl;
      if(state != ERROR)
         prevState = state;
      if(state == START) {
         if(index < argv.size() && argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            errorMessage = "No variables defined";
            state = ERROR;
         }
      }
      else if(state == VAR) {
         if(argv.size() <= index) {
            // -v but nothing after it
            errorMessage = "No variable after '-v'";
            state = ERROR;
         }
         else {
            variable = Variable::getType(argv[index]);
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-p") {
               errorMessage = "-p must be after a -d or -c";
               state = ERROR;
            }
            else {
               state = VAROPT;
            }
         }
      }
      else if(state == VAROPT) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-d") {
            state = DOWN;
            index++;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-p") {
            errorMessage = "-p must be after a -d or -c";
            state = ERROR;
         }
         else {
            // Process variable options
            vOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == NEWVAR) {
         // Check that we haven't added the variable before
         bool alreadyExists = false;
         for(int i = 0; i < variableConfigurations.size(); i++) {
            if(variableConfigurations[i].variable == variable)
               alreadyExists = true;
         }

         if(!alreadyExists) {
            dOptions.addOption("variable", Variable::getTypeName(variable));
            std::unique_ptr<Downscaler> d = Downscaler::getScheme(downscaler, variable, dOptions);
            VariableConfiguration varconf;
            varconf.variable = variable;
            varconf.downscaler = std::move(d);
            varconf.parameterFileDownscaler = std::move(parameterFileDownscaler);
            varconf.calibrators = std::move(calibrators);
            varconf.parameterFileCalibrators = std::move(parameterFileCalibrators);
            varconf.variableOptions = vOptions;
            variableConfigurations.push_back(std::move(varconf));
         }
         else {
            Util::warning("Variable '" + Variable::getTypeName(variable) + "' already read. Using first instance.");
         }

         // Reset to defaults
         vOptions.clear();
         downscaler = defaultDownscaler();
         parameterFileDownscaler = nullptr;
         dOptions.clear();
         calibrators.clear();
         parameterFileCalibrators.clear();

         if(argv.size() <= index) {
            state = END;
         }
         else {
            state = VAR;
            index++;
         }
      }
      else if(state == DOWN) {
         if(argv.size() <= index) {
            // -d but nothing after it
            errorMessage = "No downscaler after '-d'";
            state = ERROR;
         }
         else {
            downscaler = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               // Two downscalers defined for one variable
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-p") {
               state = PARDOWN;
            }
            else {
               state = DOWNOPT;
            }
         }
      }
      else if(state == DOWNOPT) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-p") {
            state = PARDOWN;
         }
         else {
            // Process downscaler options
            dOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == PARDOWN) {
         if(argv.size() <= index) {
            // -p but nothing after it
            errorMessage = "No parameter file after '-p'";
            state = ERROR;
         }
         else {
            parameterFile = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               // Two downscalers defined for one variable
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-p") {
               // Two parameter files defined for one variable
               errorMessage = "Two or more -p used for one downscaler";
               state = ERROR;
               index++;
            }
            else {
               state = PAROPTDOWN;
            }
         }
      }
      else if(state == PAROPTDOWN) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-p") {
            errorMessage = "Two or more -p used for one calibrator";
            state = ERROR;
         }
         else {
            // Process parameter file options
            pOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == CAL) {
         if(argv.size() <= index) {
            // -c but nothing after it
            state = ERROR;
            errorMessage = "No calibrator after '-c'";
         }
         else {
            calibrator = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWCAL;
            }
            else if(argv[index] == "-v") {
               state = NEWCAL;
            }
            else if(argv[index] == "-c") {
               state = NEWCAL;
            }
            else if(argv[index] == "-d") {
               state = NEWCAL;
            }
            else if(argv[index] == "-p") {
               state = PARCAL;
            }
            else {
               state = CALOPT;
            }
         }
      }
      else if(state == CALOPT) {
         if(argv.size() <= index) {
            state = NEWCAL;
         }
         else if(argv[index] ==  "-c") {
            state = NEWCAL;
         }
         else if(argv[index] == "-v") {
            state = NEWCAL;
         }
         else if(argv[index] == "-d") {
            state = NEWCAL;
         }
         else if(argv[index] == "-p") {
            state = PARCAL;
         }
         else {
            // Process calibrator options
            cOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == PARCAL) {
         if(argv.size() <= index) {
            // -p but nothing after it
            state = ERROR;
            errorMessage = "No parameter file after '-p'";
         }
         else {
            parameterFile = argv[index];
            if(argv.size() <= index) {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-c") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-d") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-p") {
               state = PARCAL;
               index++;
            }
            else {
               state = PAROPTCAL;
            }
         }
      }
      else if(state == PAROPTCAL) {
         if(argv.size() <= index) {
            state = NEWCAL;
         }
         else if(argv[index] ==  "-c") {
            state = NEWCAL;
         }
         else if(argv[index] == "-v") {
            state = NEWCAL;
         }
         else if(argv[index] == "-d") {
            state = NEWCAL;
         }
         else if(argv[index] == "-p") {
            state = PARCAL;
         }
         else {
            // Process calibrator options
            pOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == NEWCAL) {
         // We do not need to check that the same calibrator has been added for this variable
         // since this is perfectly fine (e.g. smoothing twice).

         std::unique_ptr<ParameterFile> p = nullptr;
         if(parameterFile != "") {
            p = ParameterFile::getScheme(parameterFile, pOptions);
            if(!p->isReadable()) {
               state = ERROR;
               errorMessage = "Could not open parameter file: " + pOptions.toString();
            }
         }
         if(state != ERROR) {
            cOptions.addOption("variable", Variable::getTypeName(variable));
            calibrators.push_back(Calibrator::getScheme(calibrator, cOptions));
            parameterFileCalibrators.push_back(std::move(p));

            // Reset
            calibrator = "";
            parameterFile = "";
            cOptions.clear();
            pOptions.clear();
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               state = DOWN;
               index++;
            }
            else {
               state = ERROR;
               errorMessage = "No recognized option after '-c calibrator'";
            }
         }
      }
      else if(state == END) {
         break;
      }
      else if(state == ERROR) {
         std::stringstream ss;
         ss << "Invalid command line arguments: " << errorMessage << ".";
         Util::error(ss.str());
      }
   }
}

std::string Setup::defaultDownscaler() {
   return "nearestNeighbour";
}

bool Setup::hasFile(std::string iFilename) const {
   std::map<std::string, std::unique_ptr<File>>::const_iterator it = mFileMap.find(iFilename);
   return it != mFileMap.end();
}
