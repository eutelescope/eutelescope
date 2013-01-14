/*========================================================================*/
/*          CMSPixel Decoder v2.4                                         */
/*          Author: Simon Spannagel (s.spannagel@cern.ch)                 */
/*          Created       23 feb 2012                                     */
/*          Last modified 16 aug 2012                                     */
/*========================================================================*/

#include "CMSPixelDecoder.h"

#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <map>

using namespace CMSPixel;

/*========================================================================*/
/*          CMSPixel Decoder                                              */
/*          parent class CMSPixelDecoder                                  */
/*========================================================================*/

CMSPixelDecoder::CMSPixelDecoder(const char *FileName, int *status, unsigned int rocs, int flags, unsigned int event_selection, unsigned int verbosity) 
{
    *status=0;
    CMS_stats empty = {0};
    global_statistics = empty;

    // Reading the flags:
    lazydecoding = flags & FLAG_LAZYDECODING;
    havetbm = flags & FLAG_HAVETBM;
    writeempty = flags & FLAG_EMPTYEVENTS;
    ipbus = flags & FLAG_IPBUS;
    
    // Initialize number of ROCs:
    noOfROC = rocs;
    selectevents = event_selection;
    dec_verbosity = verbosity;
    
    // Open file and check success
    log(D_LOGORRHOEA) << "Open data file..." << std::endl;
    if((mtbStream = fopen(FileName,"r")) == NULL) *status=-1;
    log(D_LOGORRHOEA) << " status: " << *status << std::endl;
}

int CMSPixelDecoder::get_event(std::vector< std::vector< CMS_event > > * event_data) {
    
    // The data vector holding one event information from the MTB file:
    std::vector< int > data;

    readmoredata:
    reject_event = false;
    event_data->clear();
    // Reset event statistics:
    pixel_addressfailed = 0;
    pixel_hits = 0;
    
    if(!chop_datastream(&data)) return 1;
    if(!process_rawdata(&data)) goto readmoredata; // return 1;
    
    unsigned int pos = 0;
    if(!check_event_sanity(&data,&pos)) goto readmoredata;

    // Init ROC id:
    unsigned int roc = 0;
    // New dummy roc vector:
    std::vector< CMS_event > evt_roc;
    // double-check if event isn't empty (maybe contains some garbage at the end, but less than a hit)
    bool have_hit = false;
    
    log(D_LOGORRHOEA) << "Looping over event data with granularity " << L_GRANULARITY << "," << L_GRANULARITY*data.size() << " iterations." << std::endl;
        
    while(pos < L_GRANULARITY*data.size()) {
        
        // Try to find a new ROC header:
        if(find_roc_header(data,&pos,roc+1)) {
            // Push back current ROC vector:
            event_data->push_back(evt_roc);
            evt_roc.clear();
            roc++;
        }
        else {
            CMS_event tmp;
            if(decode_hit(data,&pos,roc,&tmp)) {
                evt_roc.push_back(tmp);
                have_hit = true;
            }
        }
        
    }

    // Push back the last ROC vector:    
    event_data->push_back(evt_roc);
    evt_roc.clear();
    
    // If we have no hit at this point, the event is empty:
    if(!writeempty && !have_hit) {
        global_statistics.data_empty++;
        goto readmoredata;
    }
    
    // Check if the event has been rejected on the way:
    if(reject_event) {
        global_statistics.data_rejected++;
        log(D_LOGORRHOEA) << "Event rejected!" << std::endl;
        goto readmoredata;
    }
    
    // Reject events with too many ROCs at this stage, not before (they might be wrongly decoded hits...)
    if(event_data->size() > noOfROC) {
        log(D_LOGORRHOEA) << "Event contains too many ROCs!" << std::endl;
        global_statistics.data_diffrocs[event_data->size()]++;
        global_statistics.data_rejected++;
        goto readmoredata;
    }
    else if(event_data->size() < noOfROC) {
        log(D_LOGORRHOEA) << "Event contains too few ROCs!" << std::endl;
        global_statistics.data_diffrocs[event_data->size()]++;
        // selectevents == 1 means rejecting all events but those with correct No of ROC headers
        if(selectevents > 0) {
            global_statistics.data_rejected++;
            goto readmoredata;
        }
        // Push back empty ROC vectors if we allow events with less than noOfROC headers:
        else {
            global_statistics.data_fewrocs++;
            while(event_data->size() < noOfROC) event_data->push_back(evt_roc);
        }
    }
    
    if(!lazydecoding && havetbm) log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::TBM_TRAILER" << std::endl;
    
    global_statistics.data_good++;
    update_statistics();
    
    log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS end of event processing." << std::endl;
    return 0;
}

bool CMSPixelDecoder::readWord(unsigned short &word) {

    unsigned char a; // 1st byte
    unsigned char b; // 2nd byte
    
    if(feof(mtbStream) || !fread(&a,sizeof(a),1,mtbStream) || !fread(&b,sizeof(b),1,mtbStream)) {
            log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::STATUS failed to read from file. (EOF)" << std::endl;
            return false;
    }
    word = (b << 8) | a; //next word    
    log(D_LOGORRHOEA) << " " << std::hex << word << std::dec;
    return true;
}

bool CMSPixelDecoder::chop_datastream(std::vector< int > * rawdata) {

    unsigned short word;
    rawdata->clear();
    
    log(D_LOGORRHOEA) << "Chopping datastream at MTB headers..." << std::endl;
    
    if(!readWord(word)) return false;
    while (!word_is_data(word)) {
        // If header is detected read more words:
        if( word_is_header(word) ) {
            log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS header detected: " << std::hex << word << std::dec << std::endl;
            // read 3 more words after header:
            log(D_LOGORRHOEA) << "(drop:";
            for(int h = 1; h < 4; h++ ) {
                if(!readWord(word)) return false;
                log(D_LOGORRHOEA) << " " << std::hex << word << std::dec;
            }
            log(D_LOGORRHOEA) << ")" << std::endl;
        }
        else {
            global_statistics.headers_dropped++;       // Fill statistics
            log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS drop: " << std::hex << word << std::dec << std::endl;
        }
        if(!readWord(word)) return false;
    }
    
    log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS data header    : " << std::hex << word << std::dec << std::endl;
    global_statistics.data_headers++;      // Fill statistics
    // read 3 more words after header:
    for(int i = 1; i < 4; i++) if(!readWord(word)) return false;
            
    // read the data until the next MTB header arises:
    morewords:
    if(!readWord(word)) return false;
    while( !word_is_header(word) && !feof(mtbStream)){
        rawdata->push_back(word);
        if(!readWord(word)) return false;
    }

    // For IPBus readout check the next header words, too - they should be header again:
    if(!readWord(word)) return false;
    if(ipbus && !word_is_header(word)) {
        rawdata->push_back(word);
        goto morewords;
    }
    else log(D_LOGORRHOEA) << "Double-checked header." << std::endl;
    
    log(D_LOGORRHOEA) << "Raw data array size: " << rawdata->size() << ", so " << 16*rawdata->size() << " bits." << std::endl;
    
    // Rewind one word to detect the header correctly later on:
    fseek (mtbStream , -4 , SEEK_CUR);
    return true;
}

bool CMSPixelDecoder::check_event_sanity(std::vector< int > * data, unsigned int * pos) {

    log(D_LOGORRHOEA) << "Checking event sanity..." << std::endl;
    unsigned int length = L_GRANULARITY*data->size();
    
    // Checking for empty events and skip them if necessary:
    if( length <= L_EMPTYEVT ){
            global_statistics.data_empty++;        // Fill statistics
            if(!writeempty) {
                log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS event was empty, " << length << " words. Skipped." << std::endl;
                return false;
            }
            else {
                log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS event was empty, writing emtpy event." << std::endl;
                data->clear();
                return true;
            }
    }

    // Checking for data length if necessary:    
    if(!lazydecoding)
    {
        //FIXME this need to be changed in order to work with the digital ROC/TBM!
        // For digital ROC L_HIT = 24, L_ROC_HEADER = 12 -> data should be divisible by 12. (and by 6...)
        
        // (length - L_HEADER - L_TRAILER) should give a number that is divisible by 6
        // For analog ROC L_HIT = 6, L_ROC_HEADER = 3 -> data should be divisible by 6 if odd ROCs, otherwise subtract L_ROC_HEADER and try again.
        //if((((length - L_HEADER - L_TRAILER) % 6 != 0 ) && ((length - L_HEADER - L_TRAILER - L_ROC_HEADER) % 6 != 0 )) || (length - L_HEADER - L_TRAILER < 6)) {
        //FIXME look over this code, should be okay but had no time to check carefully...
        if((length - L_HEADER - L_TRAILER - noOfROC*L_ROC_HEADER) % L_HIT != 0) {
            global_statistics.data_wronglength++;      // Fill statistics
            log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::WARNING incorrect data length (" << length - L_HEADER - L_TRAILER << "). Skipped." << std::endl;
            return false;
        }
        // There might even be a huge event:
       	else if( length > L_HUGE_EVENT ) {
            global_statistics.data_huge++;     // Fill statistics
            log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::WARNING detected huge event (" << length << " words). Skipped." << std::endl;
            return false;
       	}
    }

    // If we have TBM headers check for them and remove them from data:
    // Check for missing TBM trailer if necessary:
    if(havetbm && !find_tbm_trailer(*data,length - L_TRAILER)) {
   	    if(!lazydecoding) {
            global_statistics.data_notrailer++;     // Fill statistics
            log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::WARNING: event contained no valid TBM_TRAILER. Skipped." << std::endl;
            return false;
        }
    }
    else if(L_TRAILER > 0) {
        // Delete the trailer from valid hit data:
        data->erase(data->end() - L_TRAILER,data->end());
        log(D_LOGORRHOEA) << "Deteled trailer: " << L_TRAILER << " words." << std::endl;
    }
    
    if(L_HEADER > 0) {
        // Delete the header from valid hit data:
        data->erase(data->begin(),data->begin()+L_HEADER);
        log(D_LOGORRHOEA) << "Deleted header: " << L_HEADER << " words." << std::endl;
    }

    // Maybe we deleted something, recompute length:
    length = L_GRANULARITY*data->size();

    // Find the start position with the first ROC header (maybe there are idle patterns before...)
    //FIXME the if() is an ugly hack for single analog/xdb chip readout without TBM or emulator: ROC UB level seems broken there...
    if(!havetbm && noOfROC == 1) {
        // Just set the starting position after the first ROC header:
        *pos = L_ROC_HEADER;
        log(D_LOGORRHOEA) << "Set starting position to: " << *pos << std::endl;
    }
    else {
        unsigned int count_rocs = 0;
        for(unsigned int i = 0; i < length; i++) {
            if(find_roc_header(*data, &i, 0)) {
                if(*pos == 0) {
                    log(D_LOGORRHOEA) << "Starting position after first ROC header: " << i << std::endl;
                    *pos = i;
                }
                count_rocs++;
                i--;
            }
        }
        if(count_rocs == 0) {
            log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::WARNING event contains no ROC header. Skipped." << std::endl;
            global_statistics.data_norocs++;
            return false;
        }
    /*    else if(noOfROC != count_rocs) {
            statistics.data_fewrocs++;
            statistics.data_diffrocs[count_rocs]++;
            //FIXME unsafe!
            //if(lazydecoding) {
            //  log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS ROCs: PRESET " << noOfROC << ", DATA " << count_rocs << "." << std::endl;
            //}
            // else {
            log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::WARNING ROCS: PRESET " << noOfROC << " != DATA " << count_rocs << ". Skipped." << std::endl;
            return false;
            // }
        }
        else log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS correctly detected " << count_rocs << " ROCs." << std::endl;*/
    }
    
    log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::STATUS event: " << length << " data words." << std::endl;
    for(unsigned int i = 0; i < data->size();i++) log(D_LOGORRHOEA) << data->at(i) << " ";
    log(D_LOGORRHOEA) << std::endl;
    

    return true;
}

void CMSPixelDecoder::update_statistics() {
    global_statistics.pixel_addressfailed += pixel_addressfailed;
    global_statistics.pixel_hits += pixel_hits;
}

void CMSPixelDecoder::print_statistics() {
    log(D_ERROR) << "  Processor statistics: " << std::endl;
    log(D_ERROR) << "  Detected data headers: " << global_statistics.data_headers << std::endl;
    log(D_ERROR) << "  Sane events: " << global_statistics.data_good << std::endl;
    log(D_ERROR) << "       with " << global_statistics.pixel_hits << " sane pixel hits in total." << std::endl;
    log(D_ERROR) << "       failed to convert " << global_statistics.pixel_addressfailed << " pixel addresses." << std::endl;
    log(D_ERROR) << "  Empty events: " << global_statistics.data_empty << std::endl;
    log(D_ERROR) << "  Rejected events: " << global_statistics.data_rejected << std::endl;
    log(D_ERROR) << "  Events without ROC headers: " << global_statistics.data_norocs << std::endl;
    log(D_ERROR) << "  Evaluated events with wrong # of ROC headers: " << global_statistics.data_fewrocs << std::endl;
    log(D_ERROR) << "       with ";
    for( std::map<unsigned int,int>::iterator ii=global_statistics.data_diffrocs.begin(); ii!=global_statistics.data_diffrocs.end(); ++ii)
        log(D_ERROR) << (*ii).first << " ROCs: " << (*ii).second << "x, ";
    log(D_ERROR) << std::endl;
    log(D_ERROR) << "  Variety of ROC headers: ";
    for( std::map<unsigned int,int>::iterator ii=global_statistics.data_rocheads.begin(); ii!=global_statistics.data_rocheads.end(); ++ii)
        log(D_ERROR) << std::hex << (*ii).first << " heads: " << std::dec << (*ii).second << "x, ";
    log(D_ERROR) << std::endl;
    if (!lazydecoding && havetbm) log(D_ERROR) << "  Events without TBM trailer: " << global_statistics.data_notrailer << std::endl;
    if (!lazydecoding) log(D_ERROR) << "  Events with wrong data length: " << global_statistics.data_wronglength << std::endl;
    if (!lazydecoding) log(D_ERROR) << "  Huge events: " << global_statistics.data_huge << std::endl;
    if (lazydecoding) log(D_ERROR) << "  THIS FILE WAS DECODED USING THE lazyDecoding SWITCH!" << std::endl << "  NO DETECTION OF CORRUPT EVENTS!" << std::endl;
    log(D_ERROR) << "  Event selection cut was: " << selectevents << std::endl;
    log(D_ERROR) << "  Dropped headers: " << global_statistics.headers_dropped << " (no trigger/data/reset)" << std::endl;
}

bool CMSPixelDecoder::convertDcolToCol(int dcol, int pix, int & col, int & row)
{
  log(D_LOGORRHOEA) << "converting dcol " << dcol << " pix " << pix;
  if( dcol < 0 || dcol >= ROCNUMDCOLS || pix < 2 || pix > 161) {
    row = -1;
    col = -1;
    return false;
  }

  int colEvenOdd = pix%2;  // 0 = 1st col, 1 = 2nd col of a double column

  col = dcol * 2 + colEvenOdd; // col address, starts from 0

  row = abs( int(pix/2) - ROCNUMROWS);  // row address, starts from 0

  if( col < 0 || col > ROCNUMCOLS || row < 0 || row > ROCNUMROWS ) {
    row = -1;
    col = -1;
    log(D_LOGORRHOEA) << " ...failed!" << std::endl;
    return false;
  }
  log(D_LOGORRHOEA) << " to col " << col << " row " << row << std::endl;
  return true;
}





/*========================================================================*/
/*          CMSPixel Decoder                                              */
/*          child class CMSPixelDecoderDigital                            */
/*          decoding DIGITAL chip data                                    */
/*========================================================================*/

CMSPixelDecoderDigital::CMSPixelDecoderDigital(const char *FileName, int *status, unsigned int rocs, int ana_flags, unsigned int selectevents, unsigned int verbosity) : CMSPixelDecoder(FileName, status, rocs, ana_flags, selectevents, verbosity)
{
    log(D_DEBUG) << "Preparing a digital decoder instance..." << std::endl;    
    // Loading constants:
    log(D_LOGORRHOEA) << "Loading constants..." << std::endl;
    load_constants();
}

bool CMSPixelDecoderDigital::process_rawdata(std::vector< int > * data) {

    log(D_LOGORRHOEA) << "Processing the raw data..." << std::endl;

    if(!ipbus) {
        //FIXME currently only necessary for RAL testboard environment with Altera USB readout:
        // Currently we need to dump some zero bits with the hybrid testboard setup...
        // 0000 0000 0000 0111 -> 0111
        std::vector< int > rawdata = *data;
        data->clear();

        for(unsigned int i = 0; i < rawdata.size(); i+=4) {
            unsigned short tempword = 0;
            // Take four words and combine them, each word contains 4 bits in data:
            int tail = rawdata.size() - i;
            if(tail >= 4) tail = 4;

            for(int j = 0; j < tail; j++) tempword = tempword | ((rawdata[i+j]&15) << (4*(3-j)));
            data->push_back(tempword);
        }
    }
    else {
        // IPBus data format: we need to delete some additional headers from the test board:
        
        if(data->size() > 0) {
            // remove padding to 32bit words at the end of the event by reading the data length:
            unsigned int event_length = ((data->at(1)&0xffff0000) | (data->at(0)&0x0000ffff)) - 14;
            log(D_LOGORRHOEA) << "IPBus event length: " << event_length << std::endl;
            
            // Check for stupid event sizes:
            if(event_length/2 > data->size()) return false;
            
            // cut first 8 bytes from header:
            data->erase(data->begin(),data->begin()+4);
            //  and last 14 bytes plus the padding from trailer:
            data->erase(data->begin() + (event_length/2 + event_length%2),data->end());
            
            // Swap endianness of the data:
            for(unsigned int i = 0; i < data->size(); i++) {
                unsigned int swapped = ((data->at(i)<<8)&0xff00) | ((data->at(i)>>8)&0xff);
                data->at(i) = swapped;
            }
        }
        
    }
    
    for(unsigned int i = 0; i < data->size(); i++) log(D_LOGORRHOEA) << std::hex << std::setw(4) << data->at(i) << " ";
    log(D_LOGORRHOEA) << std::dec << std::endl;
    log(D_LOGORRHOEA) << "Data array size: " << data->size() << ", so " << 16*data->size() << " bits." << std::endl;
    return true;
}

int CMSPixelDecoderDigital::get_bit(std::vector< int > data, int bit_offset) {
    unsigned int ibyte=bit_offset/L_GRANULARITY;
    int byteval;
    if (ibyte < data.size()) {
        byteval = data[ibyte];
    }
    else {
        byteval = 0;
        return 0;
    }
    // get bit, counting from most significant bit to least significant bit
    int ibit = (L_GRANULARITY-1) - (bit_offset-ibyte*L_GRANULARITY);
    int bitval = (byteval >> ibit)&1;
    return bitval;
}

int CMSPixelDecoderDigital::get_bits(std::vector< int > data, int bit_offset,int number_of_bits) {
    int value = 0;
    for (int ibit = 0; ibit < number_of_bits; ibit++) {
        value = 2*value+get_bit(data,bit_offset+ibit);
    }
    return value;
}

bool CMSPixelDecoderDigital::find_roc_header(std::vector< int > data, unsigned int * pos, unsigned int roc)
{
    // ROC header: 0111 1111 10SD, S & D are reserved status bits.
    // Possible ROC headers: 0x7f8 0x7f9 0x7fa 0x7fb
    int search_head = get_bits(data, *pos, L_ROC_HEADER);
    
    //FIXME very ugly hack to compensate for the current chip bug sending correct ROC headers: 0x7fX might be 0x3fX too...
    if(search_head == 0x7f0 || search_head == 0x3f8 || search_head == 0x3f9 || search_head == 0x3fa || search_head == 0x3fb) {
        // With selectevents == 2 only ROC headers without bit errors:
        if(selectevents > 1) reject_event = true;
        log(D_LOGORRHOEA) << "Found ROC header with bit error (" << std::hex << std::setw(4) << search_head << ") after " << std::dec << *pos << " bit." << std::endl;
        global_statistics.data_rocheads[search_head]++;
	    *pos += L_ROC_HEADER;
	    return true;
    }
    else if(search_head == 0x7f8 || search_head == 0x7f9 || search_head == 0x7fa || search_head == 0x7fb) {
        log(D_LOGORRHOEA) << "Found ROC header (" << std::hex << std::setw(4) << search_head << ") after " << std::dec << *pos << " bit." << std::endl;
        global_statistics.data_rocheads[search_head]++;
	    *pos += L_ROC_HEADER;
	    return true;
	}
	else return false;

}

bool CMSPixelDecoderDigital::find_tbm_trailer(std::vector< int > data, unsigned int pos)
{
	if (get_bits(data, pos, L_TRAILER)==0x7fa) {
	    log(D_LOGORRHOEA) << "Found TBM trailer at " << pos << "." << std::endl;
	    return true;
	}
    
    return false;
}

bool CMSPixelDecoderDigital::decode_hit(std::vector< int > data, unsigned int * pos, unsigned int roc, CMS_event * pixhit)
{
    if(L_GRANULARITY*data.size() - *pos < L_HIT) {
        log(D_LOGORRHOEA) << "Dropping " << L_GRANULARITY*data.size() - *pos << " bit at the end of the event." << std::endl;
        *pos = L_GRANULARITY*data.size();   // Set pointer to the end of data.
        return false;
    }
    int pixel_hit = get_bits(data,*pos,L_HIT);
    
    log(D_LOGORRHOEA) << "hit: " << std::hex << pixel_hit << std::dec << ": ";
        for(int i = 23; i >= 0; i--) {
            log(D_LOGORRHOEA) << ((pixel_hit>>i)&1);
            if(i==4 || i==5|| i==9|| i==12|| i==15|| i==18|| i==21) log(D_LOGORRHOEA) << ".";
        }
    log(D_LOGORRHOEA) << std::endl;

    // Double Column magic:
    //  dcol =  dcol msb        *6 + dcol lsb
    int dcol =  (pixel_hit>>21)*6 + ((pixel_hit>>18)&7);

    // Pixel ID magic:
    //  drow =  pixel msb           *36 + pixel nmsb *6 + pixel lsb
    //FIXME Beware: bug in the current ROC: row has to be inverted (~).
    int drow =  (~(pixel_hit>>15)&7)*36 + (~(pixel_hit>>12)&7)*6 + (~(pixel_hit>>9)&7);

    // pulseheight is 8 bit binary with a zero in the middle
    // to make sure we have less than eight consecutive ones.
    // pattern: XXXX 0 YYYY
    int pulseheight = ((pixel_hit>>5)&15)*16 + (pixel_hit&15);
    
    log(D_LOGORRHOEA) << "pixel dcol = " << std::setw(2) << dcol << " pixnum = " << std::setw(3) << drow << " pulseheight = " << ((pixel_hit>>5)&15) << "*16+" << (pixel_hit&15) << "= " << pulseheight << "(" << std::hex << pulseheight << ")" << std::dec << std::endl;
    
    // check the zero-bit:
    if (((pixel_hit>>4)&1) != 0) {
        pixel_addressfailed++;       // Fill statistics
        *pos += L_HIT; //jump to next hit.
        return false;
    }

    // Convert and check pixel address from double column address. Returns TRUE if address is okay:
    int col = -1;
    int row = -1;            

    if( convertDcolToCol( dcol, drow, col, row ) ) {
        pixel_hits++;        // Fill statistics
        pixhit->raw = pulseheight;
        pixhit->col = col;
        pixhit->row = row;

        log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::HIT ROC" << std::setw(2) << roc << " | pix " << pixhit->col << " " << pixhit->row << ", raw " << pixhit->raw << std::endl;
        *pos += L_HIT; //jump to next hit.
        return true;
    }

    // Else:
    pixel_addressfailed++;       // Fill statistics
    log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::WARNING failed to convert address of [" << std::hex << pixel_hit << std::dec << "]: dcol " << dcol << " drow " << drow << " (ROC" << roc << ", raw: " << pulseheight << ")" << std::endl;
    *pos += L_HIT; //jump to next hit.
    return false;
}





/*========================================================================*/
/*          CMSPixel Decoder                                              */
/*          child class CMSPixelDecoderAnalogue                           */
/*          decoding ANALOG chip data                                     */
/*========================================================================*/

CMSPixelDecoderAnalogue::CMSPixelDecoderAnalogue(const char *FileName, int *status, unsigned int rocs, int ana_flags, const char* addressFile, unsigned int selectevents, unsigned int verbosity) : CMSPixelDecoder(FileName, status, rocs, ana_flags, selectevents, verbosity)
{
    log(D_DEBUG) << "Preparing an analog decoder instance..." << std::endl;    
    // Loading constants:    
    load_constants();
    
    // Try to read the address levels file
    addressLevels = new levelset;
    if(!read_address_levels(addressFile,rocs)) *status=-2;
    else print_addresslevels();
}

bool CMSPixelDecoderAnalogue::process_rawdata(std::vector< int > * data) {
    
    log(D_LOGORRHOEA) << "Processing the raw data..." << std::endl;
    
    for(unsigned int i = 0; i < data->size(); i++) {
        data->at(i) = data->at(i) & 0x0fff;
        if( data->at(i) & 0x0800 ) data->at(i) -= 4096; //hex 800 = 2048
    }

    for(unsigned int i = 0; i < data->size();i++) log(D_LOGORRHOEA) << std::setw(4) << data->at(i) << " ";
    log(D_LOGORRHOEA) << std::endl;
    return true;
}

bool CMSPixelDecoderAnalogue::find_roc_header(std::vector< int > data, unsigned int * pos, unsigned int roc) {

    // Did we reach max number of ROCs read in from address levels file?
    if(roc >= addressLevels->ROC.size()) return false;
    
    // ROC header signature: UltraBlack, Black, lastDAC
    if(     findBin(data[*pos],2,addressLevels->ROC[roc].level) == 0 
         && findBin(data[*pos+1],2,addressLevels->ROC[roc].level) == 1 ) {
         *pos += L_ROC_HEADER;
         return true;
    }
    else return false;
}

bool CMSPixelDecoderAnalogue::find_tbm_trailer(std::vector< int > data, unsigned int pos) {
    // TBM trailer signature: UltraBlack, UltraBlack, Black, Black (+ 4 status)
    if(     findBin(data[pos],3,addressLevels->TBM.level) != 0 
         || findBin(data[pos+1],3,addressLevels->TBM.level) != 0 
         || findBin(data[pos+2],3,addressLevels->TBM.level) != 1 
         || findBin(data[pos+3],3,addressLevels->TBM.level) != 1) {
        return false;
    }
    else return true;
}

bool CMSPixelDecoderAnalogue::decode_hit(std::vector< int > data, unsigned int * pos, unsigned int roc, CMS_event * pixhit) {

    if(L_GRANULARITY*data.size() - *pos < L_HIT) {
        *pos = L_GRANULARITY*data.size();   // Set pointer to the end of data.
        return false;
    }
    
    // Find levels and translate them into DCOL and pixel addresses:
    int c1 = findBin( data[*pos+0], 5, addressLevels->address[roc].level );
    int c0 = findBin( data[*pos+1], 5, addressLevels->address[roc].level );
    int a2 = findBin( data[*pos+2], 5, addressLevels->address[roc].level );
    int a1 = findBin( data[*pos+3], 5, addressLevels->address[roc].level );
    int a0 = findBin( data[*pos+4], 5, addressLevels->address[roc].level );
    int aa = data[*pos+5];

    int dcol =  c1*6 + c0;
    int drow =  a2*36 + a1*6 + a0;

    int col = -1;
    int row = -1;

    // Convert and check pixel address from double column address. Returns TRUE if address is okay:
    if( convertDcolToCol( dcol, drow, col, row ) ) {
        pixel_hits++;        // Fill statistics
        pixhit->raw = aa;
        pixhit->col = col;
        pixhit->row = row;

        log(D_DEBUG) << "    CMSPixelDecoder::GET_EVENT::HIT ROC" << std::setw(2) << roc << " | pix " << pixhit->col << " " << pixhit->row << ", raw " << pixhit->raw << std::endl;

        *pos += L_HIT; //jump to next hit.
        return true;
    }
    else {
        pixel_addressfailed++;       // Fill statistics

        log(D_LOGORRHOEA) << "raw: " << data[*pos+0] << " " << data[*pos+1] << " " << data[*pos+2] << " " << data[*pos+3] << " " << data[*pos+4] << std::endl;
        log(D_WARNING) << "    CMSPixelDecoder::GET_EVENT::WARNING failed to convert address: dcol " << dcol << " drow " << drow << " (ROC" << roc << ", raw: " << aa << ")" << std::endl;
        *pos += L_HIT; //jump to next hit.
        return false;
    }
}


bool CMSPixelDecoderAnalogue::read_address_levels(const char* levelsFile, unsigned int rocs) {

    // Reading files with format defined by psi46expert::DecodeRawPacket::Print
    short level;
    char dummyString[100];
    char separation[100];
    levels *tempLevels = new levels();
    std::ifstream* file = new std::ifstream(levelsFile);

    if ( *file == 0 ){
        log(D_ERROR) << "    CMSPixelDecoder::READ_ADDRESS_LEVELS::ERROR cannot open the address levels file!" << std::endl;
        return false;
    }

    // Skip reading first labels:
    for (unsigned int iskip = 0; iskip < 4; iskip++ ) *file >> dummyString;

    // Read separation line:
    *file >> separation;

    // Skip reading other labels:
    for (unsigned int iskip = 0; iskip < 7; iskip++ ) *file >> dummyString;

    // Read TBM UltraBlack, black and address levels
    for (unsigned int ilevel = 0; ilevel < 8; ilevel++ ){
        *file >> level;
        addressLevels->TBM.level.push_back(level);

        if ( file->eof() || file->bad() ){
            log(D_ERROR) << "    CMSPixelDecoder::READ_ADDRESS_LEVELS::ERROR invalid format of address levels file!" << std::endl;
            return false;
        }
    }

    // Skip reading labels and separation lines:
    for (unsigned int iskip = 0; iskip < 9; iskip++ ) *file >> dummyString;

    // Read UltraBlack, black and address levels for each ROC
    // Skip reading ROC labels:
    *file >> dummyString;

    // Read file until we reach the second separation line (so EOF):    
    while((strcmp(dummyString,separation) != 0) && !file->eof() && !file->bad()) {

        // Read ROC Ultrablack and black levels:
        for (unsigned int ilevel = 0; ilevel < 3; ilevel++ ){
            if(!file->eof()) *file >> level; else goto finish;
            tempLevels->level.push_back(level);
        }
        addressLevels->ROC.push_back(*tempLevels);
        tempLevels->level.clear();
        
        // Read ROC address level encoding:
        for (unsigned int ilevel = 0; ilevel < 7; ilevel++ ){
            if(!file->eof()) *file >> level; else goto finish;
            tempLevels->level.push_back(level);
        }
        addressLevels->address.push_back(*tempLevels);
        tempLevels->level.clear();
        
        // Skip reading ROC labels:
        if(!file->eof()) *file >> dummyString; else goto finish;
         
    }

    finish:
    delete file;
    delete tempLevels;
    
    if (addressLevels->address.size() != rocs) {
            log(D_ERROR) << "    CMSPixelDecoder::READ_ADDRESS_LEVELS::ERROR No of ROCs in address levels file (" << addressLevels->address.size() << ") does not agree with GEAR geometry (" << rocs << ")!" << std::endl;
            return false;
    }
    else return true;
}

int CMSPixelDecoderAnalogue::findBin(int adc, int nlevel, std::vector< int > level ) {

    if( adc < level[0] ) return 0;
    for( int i = 0; i < nlevel; i++ ) if( adc >= level[i] && adc < level[i+1] ) return i;
    return nlevel;
}

void CMSPixelDecoderAnalogue::print_addresslevels() {

    log(D_DEBUG) << "    CMSPixelDecoder::STATUS TBM   header: ";
    std::vector<int>::iterator it;
    for ( it=addressLevels->TBM.level.begin() ; it < addressLevels->TBM.level.end(); it++ ) log(D_DEBUG) << std::setw(5) << *it << " ";
    for (unsigned int iroc = 0; iroc < addressLevels->address.size(); iroc++ ){
        log(D_DEBUG) << std::endl << "    CMSPixelDecoder::STATUS ROC" << std::setw(2) << iroc << " header: ";
        for ( it=addressLevels->ROC[iroc].level.begin() ; it < addressLevels->ROC[iroc].level.end(); it++ ) log(D_DEBUG) << std::setw(5) << *it << " ";
        log(D_DEBUG) << std::endl << "    CMSPixelDecoder::STATUS ROC" << std::setw(2) << iroc << " addr. levels: ";
        for ( it=addressLevels->address[iroc].level.begin() ; it < addressLevels->address[iroc].level.end(); it++ ) log(D_DEBUG) << std::setw(5) << *it << " ";
    }
    log(D_DEBUG) << std::endl;
}
