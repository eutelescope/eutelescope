#include "EUTelUtility.h"
#include "marlin/tinyxml.h"

using namespace std;
using namespace marlin;
using namespace eutelescope;

namespace eutelescope {

	class EUTelHistogram{
        
		private:
  	DISALLOW_COPY_AND_ASSIGN(EUTelHistogram)        // prevent users from making (default) copies of processors
      
		public:
		EUTelHistogram(std::string name, std::string);

		int book(std::string name);
		
		void setHistogramInfoName(std::string input){_histoInfoFileName = input;}

		void stringSplitting(std::string input, std::vector<std::string> &output);

		void pushStringTogether(std::string input1, std::string input2, std::string & output);


		protected:
		std::string _histoInfoFileName;
		std::string _processorName;

	};

}
