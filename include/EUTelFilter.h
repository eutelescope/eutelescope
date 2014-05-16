/* 
 *
 * \class Filter class
 * \brief Base class for object filters
 * 
 * Contains basic functionality for composition of filters
 * allowing logical operations on filters (AND, OR, NOT).
 * 
 */

#ifndef EUTELFILTER_H
#define	EUTELFILTER_H

// C++
#include <string>
#include <algorithm>
#include <vector>

namespace eutelescope {

    template < class T >
    class EUTelFilter {
        /** \typedef Member function object functor
         */
      typedef std::binder1st < std::const_mem_fun1_t< bool, EUTelFilter, T > > Binder;

    public:
        EUTelFilter() :
        _name("Filter"),
        _description("Empty filter"),
        _tester( this ) {} ;
        EUTelFilter(const EUTelFilter& orig) :
        _name( orig.GetName() ),
        _description( orig.GetDescription() ),
        _tester( &orig ) {} ;
        EUTelFilter( const std::string& name, const std::string& description ) :
        _name(name),
        _description(description),
        _tester( this ) {} ;
        
        virtual ~EUTelFilter();

        /** Test function
         * 
         * This function must return true for a given
         * object if the object satisfies selection requirements.
         * This routine must be overridden by user in order to specify the requirements.
         * 
         * @param obj Object to be tested for passage of selection conditions
         * 
         * @return Acceptance flag
         */
        virtual bool Take(const T obj) const = 0;

        virtual std::vector< T > operator()(std::vector< T >& obj) const {
	  std::vector< T > result;
	  remove_copy_if(obj.begin(), obj.end(), back_inserter(result), _tester);

	  return result;
        }

        /** Composition operator (logic <b> OR <b>)
         * 
         * Modifies selection condition such that it corresponds to
         * logical <b> OR <b> between to filters. An object will be kept if
         * it satisfy at least one filter.
         * 
         * Example of usage of this operator
         * \code{.cpp}
         *      EUTelFilter<EVENT::TrackHit*> filter1;
         *      // Set requirements for filter1
         *      EUTelFilter<EVENT::TrackHit*> filter2;
         *      // Set requirements for filter2
         *      ....
         *      filter1+=filter2;
         *      vector< EVENT::TrackHit* > allhits, selectedhits;
         *      selectedhits = filter1( allhits );   
         * \endcode
         * 
         * @param addedfilter One of the operands of the composition
         */
        virtual void operator+=(const EUTelFilter* addedfilter) {
	  _tester._logicoperand = std::bind1st(std::mem_fun(&EUTelFilter::Take), addedfilter);
            _tester._logic = Tester::OR;
        }
        
        /** Composition operator (logic <b> AND <b> )
         * 
         * Modifies selection condition such that it corresponds to
         * logical <b> AND <b> between to filters. An object will be kept if
         * it satisfy at least one filter.
         * 
         * @see operator+= for examples
         * 
         * @param addedfilter One of the operands of the composition
         */
        virtual void operator*=( const EUTelFilter* addedfilter) {
	  _tester._logicoperand = std::bind1st(std::mem_fun(&EUTelFilter::Take), addedfilter);
            _tester._logic = Tester::AND;
        }
        //	  virtual EUTelFilter operator!( const EUTelFilter& ) = 0;

        /** Filter object name setter
         * 
         * @param name
         */
        void SetName(std::string& name) {
            this->_name = name;
        }

        /** Filter object name getter
         * 
         * @return name of the filter
         */
        std::string GetName() const {
            return _name;
        }

        /** Filter description setter
         * 
         * Can be any line. Useful to elaborate description of 
         * the purpose of the filter.
         * 
         * @param description
         */
        void SetDescription(std::string description) {
            this->_description = description;
        }

        /** Filter description getter
         * 
         * @return description info
         */
        std::string GetDescription() const {
            return _description;
        }

    private:

        /** \struct Tester
         * 
         *  This is a working horse of a filter.
         */
        struct Tester {

            /** Defines a type of logical operation that is performed
             * on several filters (if specified)
             */
            enum LogicOperation {
                IDENTITY, OR, AND, NOT
            };

            /** Constructor */
            Tester( EUTelFilter* parent) : 
                _logic( IDENTITY ),
		  _logicoperand( std::bind1st(std::mem_fun(&EUTelFilter::Take), parent) ),
		  _worker( std::bind1st(std::mem_fun(&EUTelFilter::Take), parent) ) { }

            bool operator()(const T obj) {
                bool result = true;
                switch (_logic) {
                    case IDENTITY:
                        result = _worker(obj);
                        break;
                    case NOT:
                        result = !_worker(obj);
                        break;
                    case OR:
                        result = LogicOr(obj);
                        break;
                    case AND:
                        result = LogicAnd(obj);
                        break;
                    default:
                        result = true;
                }
                return !result; // reverse logic, because keep object when false in EUTelFilter
            }

            /** Performs logic <b> OR <b> between two filters */
            bool LogicOr(const T obj) {
                bool result = (_logicoperand(obj) || _worker(obj));
                return result;
            }

            /** Performs logic <b> AND <b> between two filters */
            bool LogicAnd(const T obj) {
                bool result = (_logicoperand(obj) && _worker(obj));
                return result;
            }

            /** Keeps type of logical operation */
            LogicOperation _logic;

            Binder _logicoperand;       //! logicoperand keeps second argument when combining selectors
            Binder _worker;             //! worker keeps selection operator
        };

    private:
        std::string _name;                   //! Name of the filter
        std::string _description;            //! Description info of the filter

        Tester _tester;                 // Selection working horse
    };

    template < class T >
    EUTelFilter< T >::~EUTelFilter() {
    }

} // eutelescope

#endif	/* EUTELFILTER_H */

