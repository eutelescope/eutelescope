// Version: $Id$
/*========================================================================*/
/*          CMSPixel Decoder v2.4                                         */
/*          Author: Simon Spannagel (s.spannagel@cern.ch)                 */
/*          Created       23 feb 2012                                     */
/*          Last modified 16 aug 2012                                     */
/*========================================================================*/

#ifndef CMSPIXELDECODER_H_
#define CMSPIXELDECODER_H_

#include <fstream>
#include <vector>
#include <map>
#include <iostream>

// Flags:
#define D_ERROR 1
#define D_WARNING 2
#define D_DEBUG   3
#define D_LOGORRHOEA 4

#define FLAG_EMPTYEVENTS 1
#define FLAG_LAZYDECODING 2
#define FLAG_HAVETBM 4
#define FLAG_IPBUS 8

// Sensor properties:
#define ROCNUMDCOLS 26
#define ROCNUMCOLS 52
#define ROCNUMROWS 80

namespace CMSPixel {

    // Struct for raw data readout
    typedef struct {
        int col;
        int row;
        int raw;        
    } CMS_event;

    // Struct for Decoder levels
    typedef struct {
        std::vector< int > level;
    } levels;
    
    typedef struct {
        levels TBM;
        std::vector< levels > ROC;
        std::vector< levels > address;
    } levelset;

    // Struct for readout statictics
    typedef struct {
        int data_headers;
        int data_notrailer;
        int data_wronglength;
        int data_empty;
        int data_huge;
        int data_good;
        int data_rejected;
        int data_norocs;
        int data_fewrocs;
        std::map<unsigned int,int> data_diffrocs;
        std::map<unsigned int,int> data_rocheads;
        int headers_dropped;
        int pixel_addressfailed;
        int pixel_hits;
    } CMS_stats;
    

/*========================================================================*/
/*          CMSPixel Decoder                                              */
/*          parent class CMSPixelDecoder                                  */
/*========================================================================*/
    
    class CMSPixelDecoder {
        public:
            CMSPixelDecoder(const char *FileName, int *status, unsigned int rocs, int flags, unsigned int selectevents = 0, unsigned int verbosity = 1);
	    virtual ~CMSPixelDecoder() {
	      log(D_ERROR) << "~CMSPixelDecoder done." << std::endl;
	      print_statistics();
	    }
            
            CMSPixelDecoder(const CMSPixelDecoder&); 
	    void operator=(CMSPixelDecoder const&); 

            virtual int get_event(std::vector< std::vector< CMS_event > > * event_data);
	    void print_statistics();
          	            
        protected:
            inline bool word_is_data(unsigned short word) {
                // IPBus format starts with 0xFFFFFFFF, no other headers allowed.
                if(ipbus && word == 0xFFFF) return true;
                else if(ipbus) return false;
                
                if(word == 0x8001 || word == 0x8081 || word == 0x8005) return true;
                else return false;
            };
            inline bool word_is_header(unsigned short word) {
                // IPBus format doesn't know about headers other than data headers.
                if(ipbus && word == 0xFFFF) return true;
                else if(ipbus) return false;
                
                if(word == 0x8001 || word == 0x8081 || word == 0x8005 || word == 0x8004 || word == 0x8002 || word == 0x8008 || word == 0x8010) return true;
                else return false;
            };

            // Purely virtual, to be implemented in the child classes (digital/analog):
            inline virtual void load_constants() = 0;
            virtual bool find_roc_header(std::vector< int > data, unsigned int * pos, unsigned int roc) = 0;
            virtual bool find_tbm_trailer(std::vector< int > data, unsigned int pos) = 0;
            virtual bool decode_hit(std::vector< int > data, unsigned int * pos, unsigned int roc, CMS_event * hit) = 0;
            virtual bool process_rawdata(std::vector< int > * data) = 0;
            
            // These functions are the same no matter what data format we have:
            bool check_event_sanity(std::vector< int > * data, unsigned int * pos);
          	bool convertDcolToCol(int dcol, int pix, int & col, int & row);
            
            // Common variables: length of header, trailer, data bits:
            unsigned int L_HEADER, L_TRAILER, L_EMPTYEVT, L_GRANULARITY, L_HIT, L_ROC_HEADER, L_HUGE_EVENT;
            unsigned int noOfROC, selectevents;
            bool reject_event;
            // Event statistics:
            int pixel_addressfailed, pixel_hits;

            // The flags:
            bool debug, deepdebug, lazydecoding, havetbm, writeempty, ipbus;

            // Logging:
            std::ostream& log(int debug_level) {
                static nullstream dummy;
                return debug_level > dec_verbosity ? dummy : std::cout;
            }
            CMS_stats global_statistics;
        private:
            // Data file handling:
            FILE *mtbStream;
            bool readWord(unsigned short &word);
            bool chop_datastream(std::vector< int > * data);
            void update_statistics();

            int dec_verbosity;            
            struct nullstream : std::ostream {
                nullstream() : std::ostream(0) { }
            };
    };


/*========================================================================*/
/*          CMSPixel Decoder                                              */
/*          child class CMSPixelDecoderAnalogue                           */
/*          decoding ANALOG chip data                                     */
/*========================================================================*/

    class CMSPixelDecoderAnalogue : public CMSPixelDecoder {
        public:
            CMSPixelDecoderAnalogue(const char *FileName, int *status, unsigned int rocs, int flags, const char* addressFile, unsigned int selectevents = 0, unsigned int verbosity = 1);
            virtual ~CMSPixelDecoderAnalogue(){}

            CMSPixelDecoderAnalogue(const CMSPixelDecoderAnalogue&); 
	    void operator=(CMSPixelDecoderAnalogue const&); 

        protected:
            inline void load_constants() {
                // Lenth of different tokens:
                // Analog: all values given in data words (16bit)
             	L_ROC_HEADER = 3;   // ROC header
             	L_HIT = 6;          // Hit length
             	L_GRANULARITY = 1;  // Data granularity (analog: words)

                // Check whether we should have a TBM header or not:
                if(!havetbm) {
                    L_HEADER = 1;     // FPGA header without TBM emu: 1 word;
                 	L_TRAILER = 6;    // FPGA trailer without TBM emu: 6 words;
                }
                else {
                    L_HEADER = 8;
                 	L_TRAILER = 8;
                }

                L_EMPTYEVT = L_HEADER + noOfROC*L_ROC_HEADER + L_TRAILER;                
                L_HUGE_EVENT = 2222*L_GRANULARITY; // Length of a "huge" event to be discarded
            };

            bool find_roc_header(std::vector< int > data, unsigned int * pos, unsigned int roc);
            bool find_tbm_trailer(std::vector< int > adc, unsigned int pos);
            bool decode_hit(std::vector< int > data, unsigned int * pos, unsigned int roc, CMS_event * hit);
            bool process_rawdata(std::vector< int > * data);
        private:
          	void print_addresslevels();
          	int findBin(int adc, int nlevel, std::vector< int > level);
            bool read_address_levels(const char* levelsFile, unsigned int rocs);
            levelset * addressLevels;
    };
    

/*========================================================================*/
/*          CMSPixel Decoder                                              */
/*          child class CMSPixelDecoderDigital                            */
/*          decoding DIGITAL chip data                                    */
/*========================================================================*/
    
    class CMSPixelDecoderDigital : public CMSPixelDecoder {
        public:
            CMSPixelDecoderDigital(const char *FileName, int *status, unsigned int rocs, int flags, unsigned int selectevents = 0, unsigned int verbosity = 1);
            virtual ~CMSPixelDecoderDigital(){}
        protected:
            inline void load_constants() {
                // Lenth of different tokens:
                // Digital: all values given in single bits
                L_ROC_HEADER = 12;   // ROC header
                L_HIT = 24;          // Hit length
                L_GRANULARITY = 16;  // Data granularity (digital: bits)

                // Check whether we should have a TBM header or not:
                if(!havetbm) {
                    L_HEADER = 8;
                 	L_TRAILER = 28;
                }
                else {
                    L_HEADER = 28;
                 	L_TRAILER = 28;
                }

                L_EMPTYEVT = L_HEADER + noOfROC*L_ROC_HEADER + L_TRAILER;
                L_HUGE_EVENT = 2222*L_GRANULARITY; // Length of a "huge" event to be discarded
            };
            
            bool find_roc_header(std::vector< int > data, unsigned int * pos, unsigned int roc);
            bool find_tbm_trailer(std::vector< int > adc, unsigned int pos);
            bool decode_hit(std::vector< int > data, unsigned int * pos, unsigned int roc, CMS_event * hit);
            bool process_rawdata(std::vector< int > * data);
        private:
            int get_bit(std::vector< int > data, int bit_offset);
            int get_bits(std::vector< int > data, int bit_offset,int number_of_bits);
    };
    
}
#endif /*CMSPIXELDECODER_H_*/
