#ifndef KDL_EXPRESSIONTREE_EXCEPTIONS_HPP_22309480340
#define KDL_EXPRESSIONTREE_EXCEPTIONS_HPP_22309480340


#include <stdexcept>
#include <cassert>
// CONFIGURABLE CHECKS:
//#define CHECK_CACHE
//#define STACKTRACE



#ifndef EXPRESSIONGRAPH_STACKTRACE

    #define EG_ASSERT(expr) assert(expr)
    #define EG_ASSERT_MSG(expr,msg) assert((msg,expr))


    class ExpressiongraphException : public std::logic_error 
    {
        public:
            ExpressiongraphException( const std::string& description): 
                std::logic_error(description) {}
    };

    class NotImplementedException : public ExpressiongraphException 
    {
    public:
        NotImplementedException(const char* funcname= __PRETTY_FUNCTION__): ExpressiongraphException(std::string(funcname) + " : Function not yet implemented") {
        };
    };

    class NullPointerException : public ExpressiongraphException
    {
        char msg[512];
    public:
        NullPointerException(const char* funcname= __PRETTY_FUNCTION__) : ExpressiongraphException(std::string(funcname) + " : Null pointer is given as an argument") { 
        };
    };


    class FunctionException : public ExpressiongraphException 
    {
        char msg[512];
    public:
        FunctionException(const std::string& funcname= __PRETTY_FUNCTION__) : 
            ExpressiongraphException(std::string(funcname)) { 
        };
    };

    class BodyWrongTypeException : public FunctionException 
    {
        char msg[512];
    public:
        BodyWrongTypeException(const std::string& funcname= __PRETTY_FUNCTION__) : FunctionException(std::string(funcname) + " : Type of function body does not correspond to expressiongraph function definition") { 
        };
    };

    class ArgumentWrongTypeException : public FunctionException 
    {
        char msg[512];
    public:
        ArgumentWrongTypeException(const std::string& funcname= __PRETTY_FUNCTION__) : FunctionException(std::string(funcname) + " : Type of function argument does not correspond to expressiongraph function definition") { 
        };
    };

    class ArgumentNameException : public FunctionException 
    {
        char msg[512];
    public:
        ArgumentNameException(const std::string& funcname= __PRETTY_FUNCTION__) : FunctionException(std::string(funcname) + " : Unknown name for argument") { 
        };
    };

    class ArgumentIndexException : public FunctionException 
    {
        char msg[512];
    public:
        ArgumentIndexException(const std::string& funcname= __PRETTY_FUNCTION__) : FunctionException(std::string(funcname) + " : Unknown index for argument or wrong number of arguments") { 
        };
    };

    class WrongNumberOfArgumentsException : public FunctionException 
    {
        char msg[512];
    public:
        WrongNumberOfArgumentsException(const std::string& funcname= __PRETTY_FUNCTION__) : FunctionException(std::string(funcname) + " : function has the wrong number of arguments") { 
        };
    };
#else
    #define BOOST_STACKTRACE_USE_BACKTRACE
    #define BOOST_ENABLE_ASSERT_DEBUG_HANDLER
    #include <boost/stacktrace.hpp>
    #include <boost/assert.hpp>
    #include <sstream>

    #define EG_ASSERT(expr) BOOST_ASSERT(expr)
    #define EG_ASSERT_MSG(expr,msg) BOOST_ASSERT_MSG(expr,msg)

    namespace boost {
        inline void assertion_failed_msg(char const* expr, char const* msg, char const* function, char const* /*file*/, long /*line*/) {
            std::cerr << "Expression '" << expr << "' is false in function '" << function << "': " << (msg ? msg : "<...>") << ".\n"
                << "Backtrace:\n" << boost::stacktrace::stacktrace() << '\n';
            std::abort();
        }

        inline void assertion_failed(char const* expr, char const* function, char const* file, long line) {
            ::boost::assertion_failed_msg(expr, 0 /*nullptr*/, function, file, line);
        }
    } // namespace boost

    inline std::string generate_trace() {
        using namespace std;
        ostringstream s;
        s << boost::stacktrace::stacktrace();
        return s.str();
    }


    class ExpressiongraphException : public  std::logic_error 
    {
        public:
            ExpressiongraphException(
                const std::string& description,
                const std::string& trace = generate_trace()): 
                std::logic_error("Exception thrown : "+description+"\n"+trace) {}
    };
    class NotImplementedException : public ExpressiongraphException 
    {
    public:
        NotImplementedException(): ExpressiongraphException("Function not yet implemented") {};
    };

    class NullPointerException : public ExpressiongraphException
    {
    public:
        NullPointerException(): ExpressiongraphException("Null pointer is given as an argument") {};
    };


    class FunctionException : public ExpressiongraphException
    {
    public:
        FunctionException(const std::string& description) : ExpressiongraphException(description) { 
        };
    };

    class BodyWrongTypeException : public FunctionException 
    {
    public:
        BodyWrongTypeException(): FunctionException("Type of function body does not correspond to expressiongraph function definition") {};
    };

    class ArgumentWrongTypeException : public FunctionException 
    {
    public:
        ArgumentWrongTypeException() : FunctionException("Type of function argument does not correspond to expressiongraph function definition") {};
    };

    class ArgumentNameException : public FunctionException 
    {
    public:
        ArgumentNameException() : FunctionException("Unknown name for argument") {};
    };

    class ArgumentIndexException : public FunctionException 
    {
    public:
        ArgumentIndexException() : FunctionException("Unknown index for argument or wrong number of arguments") {};
    };

    class WrongNumberOfArgumentsException : public FunctionException 
    {
    public:
        WrongNumberOfArgumentsException() : FunctionException("Function has the wrong number of arguments") {};
    };

#endif



#endif

