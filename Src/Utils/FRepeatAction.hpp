#ifndef _FPROGRESSBAR_HPP_
#define _FPROGRESSBAR_HPP_

#include <thread>
#include <chrono>

/** \brief Calls a function at time interval until stopped.
 *
 * This object will call at a defined interval any function or lambda that is
 * passed to it. The function can stop the call by returning false.
 */
class FRepeatAction {

    using function = std::function<bool()>;
    /// Stop condition
    bool _continue = true;
    
    /// The function to call
    function _f;

    /// Time between each call
    std::chrono::milliseconds _sleepTime;

    /// The thread that calls #_f
    std::thread _thread;

    /// #_thread main loop
    static void run( function f, bool& cont, std::chrono::milliseconds& sleepTime) {
        while( f() && cont ) {
            if(sleepTime == std::chrono::milliseconds{0}) {
                std::cerr << "Cannot have a null sleep time." << std::endl;
                break;
            }
            std::this_thread::sleep_for(sleepTime);
        }
    }


public:
    /// Deleted default constructor
    FRepeatAction() = delete;
    /**\brief Deleted copy constructor
     * \details A thread cannot be copied.
     */
    FRepeatAction( const FRepeatAction& other ) = delete;
    /**\brief Deleted move constructor
     * \details Calls in #run prevent moving things around.
     */
    FRepeatAction( FRepeatAction&& other ) = default;

    /**\brief Constructor
     * 
     * \param f Callable object or function to call repeatedly.
     * \param sleepTime Time interval between each call.
     */
    FRepeatAction( function f, int sleepTime ) : 
        _f(f),
        _sleepTime(sleepTime) {
        if( sleepTime == 0 ) {
            std::cerr << "Cannot have a null sleep time." << std::endl;
        } else {
            _thread = std::thread(run, _f, std::ref(_continue), std::ref(_sleepTime));
        }

    }

    /// Destructor makes a call to #stop
    ~FRepeatAction() {
        stop();
    }

    /// Stops the repeated actions and joins with #_thread.
    void stop() {
        _continue = false;
        if(_thread.joinable())
            _thread.join();
    }


};




#endif
