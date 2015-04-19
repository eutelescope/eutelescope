/*
 *This is the template functions for EUTelState class.
 *The idea is that we will pass any time of precision variable and type and store it in a double lcio object. 
 *Must be template so we can interact with more generic input.
 *       
 */
#ifndef EUTELSTATE_TCC  
#define EUTELSTATE_TCC



template<class number>
void EUTelState::setPositionLocal(number position[]){
        const float pos[3] = {position[0],position[1],position[2]};
        setReferencePoint(pos);
}
#endif
