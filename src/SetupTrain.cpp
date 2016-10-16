#include "SetupTrain.h"
#include "File/File.h"
#include "Calibrator/Calibrator.h"
#include "Downscaler/Downscaler.h"

SetupTrain::SetupTrain(const std::vector<std::string>& argv) {

   if(argv.size() <= 0) {
      Util::error("Missing training data file");
   }

   // Set initial non-working values
   variable = Variable::None;

   // Implement a finite state machine
   enum State {START = 0, OUTPUT = 1, OUTPUTOPT = 2, METHOD = 3, METHODOPT = 10, VAR = 20, END = 90, ERROR = 100};
   State state = START;
   State prevState = START;

   downscaler = Downscaler::getScheme("nearestNeighbour", Variable::T, Options());

   // Process obs/fcst filenames and options
   Options obsOptions, fcstOptions;
   std::string obsFilename, fcstFilename;
   int index = 0;
   while(index < argv.size()) {
      std::string arg = argv[index];
      if(obsFilename == "") {
         obsFilename = arg;
      }
      else if(fcstFilename == "") {
         if(Util::hasChar(arg, '=')) {
            obsOptions.addOptions(arg);
         }
         else {
            fcstFilename = arg;
         }
      }
      else {
         if(Util::hasChar(arg, '=')) {
            fcstOptions.addOptions(arg);
         }
         else {
            break;
         }
      }
      index++;
   }
   std::vector<std::string> obsFilenames = Util::glob(obsFilename);
   std::vector<std::string> fcstFilenames = Util::glob(fcstFilename);

   if(obsFilenames.size() == 0 || fcstFilenames.size() == 0) {
      Util::error("No valid obs or fcst files");
   }
   for(int i = 0; i < obsFilenames.size(); i++) {
      std::unique_ptr<File> file = File::getScheme(obsFilenames[i], obsOptions, true);
      if(file == NULL) {
         Util::error("Could not open '" + obsFilenames[i] + "'");
      }
      observations.push_back(std::move(file));
   }
   for(int i = 0; i < fcstFilenames.size(); i++) {
      std::unique_ptr<File> file = File::getScheme(fcstFilenames[i], fcstOptions, true);
      if(file == NULL) {
         Util::error("Could not open '" + fcstFilenames[i] + "'");
      }
      forecasts.push_back(std::move(file));
   }

   Options oOptions;
   Options mOptions;
   std::string outputType = "";
   std::string methodType = "";
   std::string errorMessage = "";
   while(true) {
      // std::cout << state << std::endl;
      if(state != ERROR)
         prevState = state;
      if(state == START) {
         if(argv.size() <= index) {
            state = ERROR;
            errorMessage = "'-p', '-v', and '-c' required";
         }
         else if(argv[index] == "-p") {
            state = OUTPUT;
            index++;
         }
         else if(argv[index] == "-c") {
            state = METHOD;
            index++;
         }
         else if(argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            state = ERROR;
            errorMessage = "'" + argv[index] + "' unrecognized";
         }
      }
      else if(state == OUTPUT) {
         if(argv.size() <= index) {
            // -p but nothing after it
            errorMessage = "Nothing after '-p'";
            state = ERROR;
         }
         else {
            if(outputType != "") {
               state = ERROR;
               errorMessage = "Duplicate '-p'";
            }
            else {
               outputType = argv[index];
               index++;
               state = OUTPUTOPT;
            }
         }
      }
      else if(state == OUTPUTOPT) {
         if(argv.size() <= index) {
            state = END;
         }
         else if(argv[index] == "-c") {
            state = METHOD;
            index++;
         }
         else if(argv[index] == "-p") {
            state = ERROR;
            errorMessage = "Duplicate '-p'";
         }
         else if(argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            // Process output options
            oOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == METHOD) {
         if(argv.size() <= index) {
            // -c but nothing after it
            errorMessage = "Nothing after '-c'";
            state = ERROR;
         }
         else {
            if(methodType != "") {
               state = ERROR;
               errorMessage = "Duplicate '-c'";
            }
            else {
               methodType = argv[index];
               index++;
               state = METHODOPT;
            }
         }
      }
      else if(state == METHODOPT) {
         if(argv.size() <= index) {
            state = END;
         }
         else if(argv[index] == "-p") {
            state = OUTPUT;
            index++;
         }
         else if(argv[index] == "-c") {
            state = ERROR;
            errorMessage = "Duplicate '-c'";
         }
         else if(argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            // Process method options
            mOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == VAR) {
         if(argv.size() <= index) {
            // -v but nothing after it
            errorMessage = "Nothing after '-v'";
            state = ERROR;
         }
         else {
            if(variable != Variable::None) {
               state = ERROR;
               errorMessage = "Duplicate '-v'";
            }
            else {
               variable = Variable::getType(argv[index]);
               index++;
               if(argv.size() <= index) {
                  state = END;
               }
               else if(argv[index] == "-c") {
                  state = METHOD;
                  index++;
               }
               else if(argv[index] == "-p") {
                  state = OUTPUT;
                  index++;
               }
               else if(argv[index] == "-v") {
                  state = ERROR;
                  errorMessage = "Duplicate '-v'";
               }
               else {
                  state = ERROR;
                  errorMessage = "Junk after -v <variable>";
               }
            }
         }
      }
      else if(state == END) {
         mOptions.addOption("variable", Variable::getTypeName(variable));
         method = Calibrator::getScheme(methodType, mOptions);
         output = ParameterFile::getScheme(outputType, oOptions, true);
         break;
      }
      else if(state == ERROR) {
         std::stringstream ss;
         ss << "Could not understand command line arguments: " << errorMessage << ".";
         Util::error(ss.str());
      }
      else {
         abort();
      }
   }
   if(method) {
      std::stringstream ss;
      ss << "Could not understand command line arguments: Missing -c.";
      Util::error(ss.str());
   }
   if(!output) {
      std::stringstream ss;
      ss << "Could not understand command line arguments: Missing -p.";
      Util::error(ss.str());
   }
}
SetupTrain::~SetupTrain() {
   delete downscaler;
}
