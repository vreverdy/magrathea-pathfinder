/* ********************************** TIMER ********************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          Timer
// DESCRIPTION :    A timer to manage time measurements and benchmarks
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           timer.h
/// \brief          A timer to manage time measurements and benchmarks
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <functional>
#include <ratio>
#include <chrono>
#include <ctime>
#include <cmath>
#include <utility>
// Include libs
// Include project
#include "timer.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// A timer to manage time measurements and benchmarks
/// \brief          A timer to manage time measurements and benchmarks.
/// \details        Provides a wrapper of <tt>std::chrono</tt> for an easy use 
///                 and basic operations needed by time execution management. It
///                 has two internal times points : one to mark the beginning of
///                 a measurement, and one to mark the end of the current 
///                 measurement. It has also a reference time point to evaluate
///                 real elasped time.
/// \tparam         Type Duration representation type.
/// \tparam         Period Standard ratio representing the tick period.
/// \tparam         Clock Internal clock type.
template <typename Type = double, class Period = std::chrono::seconds::period, class Clock = std::chrono::steady_clock>
class Timer final
{
    // Setup
    static_assert(std::is_arithmetic<Type>::value, "ERROR = Timer : non arithmetic type");
    static_assert(Period::den != 0, "ERROR = Timer : non compatible period type");
    static_assert(Clock::period::den != 0, "ERROR = Timer : non compatible clock type");
    
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        template <class Duration = Type, class ReferenceTimePoint = std::true_type, class BeginningTimePoint = std::true_type, class EndingTimePoint = std::true_type> inline Timer(const bool running0 = false, const Duration& record0 = Duration(), const ReferenceTimePoint& reference0 = ReferenceTimePoint(), const BeginningTimePoint& beginning0 = BeginningTimePoint(), const EndingTimePoint& ending0 = EndingTimePoint());
        template <typename OtherType, class OtherPeriod> inline Timer(const Timer<OtherType, OtherPeriod, Clock>& source);
    //@}

    // Operators
    /// \name           Operators
    //@{
    public:
        template <typename OtherType, class OtherPeriod> inline Timer<Type, Period, Clock>& operator=(const Timer<OtherType, OtherPeriod, Clock>& rhs);
        inline Type operator()() const;
    //@}

    // Assignment
    /// \name           Assignment
    //@{
    public:
        template <class Duration = Type, class ReferenceTimePoint = std::true_type, class BeginningTimePoint = std::true_type, class EndingTimePoint = std::true_type> inline Timer<Type, Period, Clock>& assign(const bool running0 = false, const Duration& record0 = Duration(), const ReferenceTimePoint& reference0 = ReferenceTimePoint(), const BeginningTimePoint& beginning0 = BeginningTimePoint(), const EndingTimePoint& ending0 = EndingTimePoint());
        template <typename OtherType, class OtherPeriod> inline Timer<Type, Period, Clock>& assign(const Timer<OtherType, OtherPeriod, Clock>& source);
    //@}
    
    // Management
    /// \name           Management
    //@{
    public:
        inline Timer<Type, Period, Clock> copy() const;
        template <typename OtherType = Type, class OtherPeriod = Period> inline Timer<OtherType, OtherPeriod, Clock> cast() const;
    //@}
    
    // Getters
    /// \name           Getters
    //@{
    public:
        inline const bool& running() const;
        inline const typename Clock::duration& record() const;
        inline const typename Clock::time_point& reference() const;
        inline const typename Clock::time_point& beginning() const;
        inline const typename Clock::time_point& ending() const;
    //@}

    // Actions
    /// \name           Actions
    //@{
    public:
        inline std::chrono::duration<Type, Period> reset();
        inline std::chrono::duration<Type, Period> start();
        inline std::chrono::duration<Type, Period> stop();
    //@}

    // Measurement
    /// \name           Measurement
    //@{
    public:
        inline std::chrono::duration<Type, Period> current() const;
        inline std::chrono::duration<Type, Period> total() const;
        inline std::chrono::duration<Type, Period> real() const;
    //@}
    
    // Stream
    /// \name           Stream
    //@{
    public:
        template <typename SelfType, class SelfPeriod, class SelfClock> friend std::ostream& operator<<(std::ostream& lhs, const Timer<SelfType, SelfPeriod, SelfClock>& rhs);
    //@}
    
    // Utilities
    /// \name           Utilities
    //@{
    public:
        template <class Duration = Type, typename Counter = unsigned long long int> static inline std::chrono::duration<Type, Period> wait(const Duration& delay = Duration(1), Counter&& counter = Counter());
        template <typename Counter, class Function, class... Args> static inline std::chrono::duration<Type, Period> benchmark(const Counter& counter, Function&& f, Args&&... args);
    //@}
    
    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}

    // Data members
    /// \name           Data members
    //@{
    protected:
        bool _running;                                                          ///< Flag to indicate whether a measurement is running or not.
        typename Clock::duration _record;                                       ///< Internal backup of duration.
        typename Clock::time_point _reference;                                  ///< Reference time point for measurements.
        typename Clock::time_point _beginning;                                  ///< Beginning time point for measurements.
        typename Clock::time_point _ending;                                     ///< Ending time point for measurements.
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Implicit contents constructor
/// \brief          Implicit contents constructor.
/// \details        Provides a construction from every single parameter of the
///                 timer : the current record of duration, the time of 
///                 reference, the time of the last start, the time of the last 
///                 stop, and whether the timer is running. If no reference is 
///                 provided, the current time is set. If no beginning time is
///                 provided, it is set to the reference time. If no ending time
///                 is provided, it is set to the beginning.
/// \tparam         Duration (Standard duration type or arithmetic type.)
/// \tparam         ReferenceTimePoint (Time point type of the reference.)
/// \tparam         BeginningTimePoint (Time point type of the beginning time.)
/// \tparam         EndingTimePoint (Time point type of the ending time.)
/// \param[in]      running0 Input value of whether the timer is running.
/// \param[in]      record0 Input value of the backup of duration.
/// \param[in]      reference0 Input value of the reference time.
/// \param[in]      beginning0 Input value of the beginning time.
/// \param[in]      ending0 Input value of the ending time.
template <typename Type, class Period, class Clock>
template <class Duration, class ReferenceTimePoint, class BeginningTimePoint, class EndingTimePoint> 
inline Timer<Type, Period, Clock>::Timer(const bool running0, const Duration& record0, const ReferenceTimePoint& reference0, const BeginningTimePoint& beginning0, const EndingTimePoint& ending0)
: _running(running0)
, _record(std::chrono::duration_cast<typename Clock::duration>(typename std::conditional<std::is_arithmetic<Duration>::value, std::chrono::duration<Type, Period>, Duration>::type(record0)))
, _reference((std::is_same<ReferenceTimePoint, std::true_type>::value) ? (Clock::now()) : (std::chrono::time_point_cast<typename Clock::duration>(typename std::conditional<std::is_same<ReferenceTimePoint, std::true_type>::value, typename Clock::time_point, ReferenceTimePoint>::type(typename std::conditional<std::is_same<ReferenceTimePoint, std::true_type>::value, typename Clock::duration, ReferenceTimePoint>::type(reference0)))))
, _beginning((std::is_same<BeginningTimePoint, std::true_type>::value) ? (_reference) : (std::chrono::time_point_cast<typename Clock::duration>(typename std::conditional<std::is_same<BeginningTimePoint, std::true_type>::value, typename Clock::time_point, BeginningTimePoint>::type(typename std::conditional<std::is_same<BeginningTimePoint, std::true_type>::value, typename Clock::duration, BeginningTimePoint>::type(beginning0)))))
, _ending((std::is_same<EndingTimePoint, std::true_type>::value) ? (_beginning) : (std::chrono::time_point_cast<typename Clock::duration>(typename std::conditional<std::is_same<EndingTimePoint, std::true_type>::value, typename Clock::time_point, EndingTimePoint>::type(typename std::conditional<std::is_same<EndingTimePoint, std::true_type>::value, typename Clock::duration, EndingTimePoint>::type(ending0)))))
{
    ;
}

// Implicit conversion constructor
/// \brief          Implicit conversion constructor.
/// \details        Provides a construction from a timer of another type.
/// \tparam         OtherType (Other duration representation type.)
/// \tparam         OtherPeriod (Other standard ratio representing the tick
///                 period.)
/// \param[in]      source Source of the copy.
template <typename Type, class Period, class Clock>
template <typename OtherType, class OtherPeriod> 
inline Timer<Type, Period, Clock>::Timer(const Timer<OtherType, OtherPeriod, Clock>& source)
: _running(source.running())
, _record(source.record())
, _reference(source.reference())
, _beginning(source.beginning())
, _ending(source.ending())
{
    ;
}
// -------------------------------------------------------------------------- //



// -------------------------------- OPERATORS ------------------------------- //
// Conversion assignment operator
/// \brief          Conversion assignment operator.
/// \details        Assigns the contents from a timer of another type.
/// \tparam         OtherType (Other duration representation type.)
/// \tparam         OtherPeriod (Other standard ratio representing the tick
///                 period.)
/// \param[in]      rhs Right-hand side.
/// \return         Self reference.
template <typename Type, class Period, class Clock>
template <typename OtherType, class OtherPeriod> 
inline Timer<Type, Period, Clock>& Timer<Type, Period, Clock>::operator=(const Timer<OtherType, OtherPeriod, Clock>& rhs)
{
    _running = rhs.running();
    _record = rhs.record();
    _reference = rhs.reference();
    _beginning = rhs.beginning();
    _ending = rhs.ending();
    return *this;
}

// Total duration extraction operator
/// \brief          Total duration extraction operator.
/// \details        Computes the total elapsed duration between all starts and
///                 stops since the last reset and convert it to an arithmetic
///                 type.
/// \return         Total duration.
template <typename Type, class Period, class Clock>
inline Type Timer<Type, Period, Clock>::operator()() const
{
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >((_running) ? (_record+(Clock::now()-_beginning)) : (_record)).count();
}
// -------------------------------------------------------------------------- //



// ------------------------------- ASSIGNMENT ------------------------------- //
// Contents assignment
/// \brief          Contents assignment.
/// \details        Assigns contents from every single parameter of the timer : 
///                 the current record of duration, the time of reference, the 
///                 time of the last start, the time of the last stop, and whether 
///                 the timer is running. If no reference is provided, the current 
///                 time is set. If no beginning time is provided, it is set to 
///                 the reference time. If no ending time is provided, it is set 
///                 to the beginning.
/// \tparam         Duration (Standard duration type or arithmetic type.)
/// \tparam         ReferenceTimePoint (Time point type of the reference.)
/// \tparam         BeginningTimePoint (Time point type of the beginning time.)
/// \tparam         EndingTimePoint (Time point type of the ending time.)
/// \param[in]      running0 Input value of whether the timer is running.
/// \param[in]      record0 Input value of the backup of duration.
/// \param[in]      reference0 Input value of the reference time.
/// \param[in]      beginning0 Input value of the beginning time.
/// \param[in]      ending0 Input value of the ending time.
/// \return         Self reference.
template <typename Type, class Period, class Clock>
template <class Duration, class ReferenceTimePoint, class BeginningTimePoint, class EndingTimePoint> 
inline Timer<Type, Period, Clock>& Timer<Type, Period, Clock>::assign(const bool running0, const Duration& record0, const ReferenceTimePoint& reference0, const BeginningTimePoint& beginning0, const EndingTimePoint& ending0)
{
    _running = running0;
    _record = std::chrono::duration_cast<typename Clock::duration>(typename std::conditional<std::is_arithmetic<Duration>::value, std::chrono::duration<Type, Period>, Duration>::type(record0));
    _reference = (std::is_same<ReferenceTimePoint, std::true_type>::value) ? (Clock::now()) : (std::chrono::time_point_cast<typename Clock::duration>(typename std::conditional<std::is_same<ReferenceTimePoint, std::true_type>::value, typename Clock::time_point, ReferenceTimePoint>::type(typename std::conditional<std::is_same<ReferenceTimePoint, std::true_type>::value, typename Clock::duration, ReferenceTimePoint>::type(reference0))));
    _beginning = (std::is_same<BeginningTimePoint, std::true_type>::value) ? (_reference) : (std::chrono::time_point_cast<typename Clock::duration>(typename std::conditional<std::is_same<BeginningTimePoint, std::true_type>::value, typename Clock::time_point, BeginningTimePoint>::type(typename std::conditional<std::is_same<BeginningTimePoint, std::true_type>::value, typename Clock::duration, BeginningTimePoint>::type(beginning0))));
    _ending = (std::is_same<EndingTimePoint, std::true_type>::value) ? (_beginning) : (std::chrono::time_point_cast<typename Clock::duration>(typename std::conditional<std::is_same<EndingTimePoint, std::true_type>::value, typename Clock::time_point, EndingTimePoint>::type(typename std::conditional<std::is_same<EndingTimePoint, std::true_type>::value, typename Clock::duration, EndingTimePoint>::type(ending0))));
    return *this;
}

// Conversion assignment
/// \brief          Conversion assignment.
/// \details        Assigns the contents from a timer of another type.
/// \tparam         OtherType (Other duration representation type.)
/// \tparam         OtherPeriod (Other standard ratio representing the tick
///                 period.)
/// \param[in]      source Source of the copy.
/// \return         Self reference.
template <typename Type, class Period, class Clock>
template <typename OtherType, class OtherPeriod> 
inline Timer<Type, Period, Clock>& Timer<Type, Period, Clock>::assign(const Timer<OtherType, OtherPeriod, Clock>& source)
{
    _running = source.running();
    _record = source.record();
    _reference = source.reference();
    _beginning = source.beginning();
    _ending = source.ending();
    return *this;
}
// -------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Copy
/// \brief          Copy.
/// \details        Returns a copy of the timer.
/// \return         Copy.
template <typename Type, class Period, class Clock>
inline Timer<Type, Period, Clock> Timer<Type, Period, Clock>::copy() const
{
    return *this;
}

// Cast
/// \brief          Cast.
/// \details        Casts the timer to another timer type with another period.
/// \tparam         OtherType Other duration representation type.
/// \tparam         OtherPeriod Other standard ratio representing the tick 
///                 period.
/// \return         Copy.
template <typename Type, class Period, class Clock>
template <typename OtherType, class OtherPeriod> 
inline Timer<OtherType, OtherPeriod, Clock> Timer<Type, Period, Clock>::cast() const
{
    return Timer<OtherType, OtherPeriod, Clock>(*this);
}
// -------------------------------------------------------------------------- //



// --------------------------------- GETTERS -------------------------------- //
// Get whether the timer is running
/// \brief          Get whether the timer is running.
/// \details        Returns true if the timer is running, false if it was 
///                 stopped.
/// \return         Running status.
template <typename Type, class Period, class Clock>
inline const bool& Timer<Type, Period, Clock>::running() const
{
    return _running;
}

// Get the saved duration
/// \brief          Get the saved duration.
/// \details        Returns the value of the saved duration which is the total 
///                 duration saved during the last stop.
/// \return         Duration.
template <typename Type, class Period, class Clock>
inline const typename Clock::duration& Timer<Type, Period, Clock>::record() const
{
    return _record;
}

// Get the reference time point
/// \brief          Get the reference time point.
/// \details        Returns the reference time point, which is generally the 
///                 time point of the last reset.
/// \return         Reference time point.
template <typename Type, class Period, class Clock>
inline const typename Clock::time_point& Timer<Type, Period, Clock>::reference() const
{
    return _reference;
}

// Get the beginning time point
/// \brief          Get the beginning time point.
/// \details        Returns the beginning time point, which is generally the 
///                 time point of the last start.
/// \return         Beginning time point.
template <typename Type, class Period, class Clock>
inline const typename Clock::time_point& Timer<Type, Period, Clock>::beginning() const
{
    return _beginning;
}

// Get the ending time point
/// \brief          Get the ending time point.
/// \details        Returns the ending time point, which is generally the 
///                 time point of the last stop. If the timer is running, then
///                 the ending time point is equal to the last start time point.
/// \return         Ending time point.
template <typename Type, class Period, class Clock>
inline const typename Clock::time_point& Timer<Type, Period, Clock>::ending() const
{
    return _ending;
}
// -------------------------------------------------------------------------- //



// --------------------------------- ACTIONS -------------------------------- //
// Reset timer
/// \brief          Reset timer.
/// \details        Resets the timer : all time points are set to the current
///                 time, the duration is set to zero, and the timer is set off.
/// \return         Current duration, which is equal to zero.
template <typename Type, class Period, class Clock>
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::reset()
{
    _running = false;
    _record = Clock::duration::zero();
    _reference = Clock::now();
    _beginning = _reference;
    _ending = _reference;
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >(_record);
}

// Start timer
/// \brief          Start timer.
/// \details        Starts the timer for a new measurement. If the timer is 
///                 already running, the previous state is erased.
/// \return         Current duration, which is equal to zero.
template <typename Type, class Period, class Clock>
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::start()
{
    _beginning = Clock::now();
    _ending = _beginning;
    _running = true;
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >(_ending-_beginning);
}

// Stop timer
/// \brief          Stop timer.
/// \details        Stops the timer, adds the duration to the total one, and 
///                 returns the time since the previous start. If the timer is
///                 already not running, nothing is done.
/// \return         Current duration.
template <typename Type, class Period, class Clock>
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::stop()
{
    if (_running) {
        _ending = Clock::now();
        _record += _ending-_beginning;
        _running = false;
    }
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >(_ending-_beginning);
}
// -------------------------------------------------------------------------- //



// ------------------------------ MEASUREMENT ------------------------------- //
// Duration since last start
/// \brief          Duration since last start.
/// \details        Computes the duration since last start : if the timer is 
///                 still running, it computes the difference between the call
///                 time and the last start time and if the timer is not running
///                 it returns the difference between the last start and the
///                 last stop.
/// \return         Current duration.
template <typename Type, class Period, class Clock>
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::current() const
{
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >((_running) ? (Clock::now()-_beginning) : (_ending-_beginning));
}

// Total duration
/// \brief          Total duration.
/// \details        Computes the total elapsed duration between all starts and
///                 stops since the last reset.
/// \return         Total duration.
template <typename Type, class Period, class Clock>
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::total() const
{
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >((_running) ? (_record+(Clock::now()-_beginning)) : (_record));
}

// Real duration
/// \brief          Real duration.
/// \details        Computes the real duration since the last reset without 
///                 considering any start and stop.
/// \return         Real duration.
template <typename Type, class Period, class Clock>
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::real() const
{
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >(Clock::now()-_reference);
}
// -------------------------------------------------------------------------- //



// --------------------------------- STREAM --------------------------------- //
// Output stream operator
/// \brief          Output stream operator.
/// \details        Prints out the total duration.
/// \tparam         SelfType (Duration representation type.)
/// \tparam         SelfPeriod (Standard ratio representing the tick period.)
/// \tparam         SelfClock (Internal clock type.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side timer.
/// \return         Output stream.
template <typename SelfType, class SelfPeriod, class SelfClock>
std::ostream& operator<<(std::ostream& lhs, const Timer<SelfType, SelfPeriod, SelfClock>& rhs)
{
    return lhs<<std::chrono::duration_cast<std::chrono::duration<SelfType, SelfPeriod> >((rhs._running) ? (rhs._record+(SelfClock::now()-rhs._beginning)) : (rhs._record)).count();
}
// -------------------------------------------------------------------------- //



// ------------------------------- UTILITIES -------------------------------- //
// Wait a certain time
/// \brief          Wait a certain time.
/// \details        Loops over time in order to delay some operation. If a 
///                 counter is passed, it is incremented at each loop. The loop
///                 ends when the elapsed time is greater or equal to the 
///                 specified delay. This delay can be a number or a standard 
///                 duration. Internally, it is casted to the timer type, and 
///                 then to a high precision clock. At the end of the function
///                 the total duration is casted back from the high precision 
///                 clock to the timer duration type.
/// \tparam         Duration (Standard duration type or arithmetic type.)
/// \tparam         Counter (Type that can be incremented.)
/// \param[in]      delay Duration to wait.
/// \param[in,out]  counter Incremented counter at each loop step.
/// \return         The real elapsed time in the function.
template <typename Type, class Period, class Clock>
template <class Duration, typename Counter> 
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::wait(const Duration& delay, Counter&& counter)
{
    const typename Clock::duration limit = std::chrono::duration_cast<typename Clock::duration>(typename std::conditional<std::is_arithmetic<Duration>::value, std::chrono::duration<Type, Period>, Duration>::type(delay));
    const typename Clock::time_point marker = Clock::now();
    while (Clock::now()-marker < limit) {
        ++counter;
    }
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >(Clock::now()-marker);
}

// Benckmark a function
/// \brief          Benchmark a function.
/// \details        Executes the provided function in a loop and computes the
///                 total time needed to run it. The call uses a volatile 
///                 temporary to prevent null statement optimization, but some
///                 compilers may manage to optimize that. The results are not
///                 guaranteed to be exact and should be checked with a real
///                 benchmarking suite. The returned time correspond to the time
///                 of function execution plus the time of the copy.
/// \tparam         Counter (Type that can be incremented.)
/// \tparam         Function (Function type : <tt>Something(Args...)</tt>.)
/// \tparam         Args (Arguments types.) 
/// \param[in]      counter Number of loops to do.
/// \param[in]      f Function object <tt>Something(Args...)</tt>.
/// \param[in]      args Arguments of the function.
/// \return         Total duration of loop execution.
template <typename Type, class Period, class Clock>
template <typename Counter, class Function, class... Args>
inline std::chrono::duration<Type, Period> Timer<Type, Period, Clock>::benchmark(const Counter& counter, Function&& f, Args&&... args)
{
    volatile decltype(f(args...)) temporary = decltype(f(args...))();
    const typename Clock::time_point marker = Clock::now();
    for (Counter i = Counter(); i < counter; ++i) {
        temporary = f(std::forward<Args>(args)...);
    }
    *reinterpret_cast<volatile unsigned char*>(&temporary) = *reinterpret_cast<volatile unsigned char*>(&temporary) + static_cast<unsigned char>(0);
    return std::chrono::duration_cast<std::chrono::duration<Type, Period> >(Clock::now()-marker);
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Timer.
/// \return         0 if no error.
template <typename Type, class Period, class Clock>
int Timer<Type, Period, Clock>::example()
{
    // Initialize
    std::cout<<"BEGIN = Timer::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;

    // Construction
    Timer<double> ts(true);
    Timer<double, std::chrono::milliseconds::period> tms;
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"Timer<double>() : "                                                                 <<Timer<double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"Timer<double>(false, 10, Clock::now()) : "                                          <<Timer<double>(false, 10, Clock::now())<<std::endl;
    std::cout<<std::setw(width*2)<<"Timer<double>(true, 10, Clock::now()) : "                                           <<Timer<double>(true, 10, Clock::now())<<std::endl;
    std::cout<<std::setw(width*2)<<"Timer<long long int, std::chrono::milliseconds::period>() : "                       <<Timer<long long int, std::chrono::milliseconds::period>()<<std::endl;
    std::cout<<std::setw(width*2)<<"Timer<long long int, std::chrono::milliseconds::period>(ts) : "                     <<Timer<long long int, std::chrono::milliseconds::period>(ts)<<std::endl;

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"tms = ts : "                                  <<(tms = ts)<<std::endl;
    std::cout<<std::setw(width)<<"ts()"                                         <<ts()<<std::endl;
    
    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"tms.assign() : "                              <<tms.assign()<<std::endl;
    std::cout<<std::setw(width)<<"tms.assign(false, 10) : "                     <<tms.assign(false, 10)<<std::endl;
    std::cout<<std::setw(width)<<"tms.assign(true, 10) : "                      <<tms.assign(true, 10)<<std::endl;
    std::cout<<std::setw(width)<<"tms.assign(ts) : "                            <<tms.assign(ts)<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"tms.copy() : "                                <<tms.copy()<<std::endl;
    std::cout<<std::setw(width)<<"tms.cast() : "                                <<tms.cast()<<std::endl;
    std::cout<<std::setw(width)<<"tms.cast<int, std::ratio<1> >() : "           <<tms.cast<int, std::ratio<1> >()<<std::endl;

    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Getters : "                                                                         <<std::endl;
    std::cout<<std::setw(width*2)<<"ts.running() : "                                                                    <<ts.running()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.record().count() : "                                                             <<ts.record().count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.reference().time_since_epoch().count() : "                                       <<ts.reference().time_since_epoch().count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.beginning().time_since_epoch().count() : "                                       <<ts.beginning().time_since_epoch().count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.ending().time_since_epoch().count() : "                                          <<ts.ending().time_since_epoch().count()<<std::endl;

    // Actions
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Actions : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"tms.reset().count() : "                       <<tms.reset().count()<<std::endl;
    std::cout<<std::setw(width)<<"tms.start().count() : "                       <<tms.start().count()<<std::endl;
    std::cout<<std::setw(width)<<"tms.stop().count() : "                        <<tms.stop().count()<<std::endl;
    std::cout<<std::setw(width)<<"tms.start().count() : "                       <<tms.start().count()<<std::endl;
    std::cout<<std::setw(width)<<"tms.stop().count() : "                        <<tms.stop().count()<<std::endl;
    std::cout<<std::setw(width)<<"tms.stop().count() : "                        <<tms.stop().count()<<std::endl;

    // Measurement
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Measurement : "                               <<std::endl;
    std::cout<<std::setw(width)<<"tms.current().count() : "                     <<tms.current().count()<<std::endl;
    std::cout<<std::setw(width)<<"tms.total().count() : "                       <<tms.total().count()<<std::endl;
    std::cout<<std::setw(width)<<"tms.real().count() : "                        <<tms.real().count()<<std::endl;

    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, ts) : "                 <<ts<<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, tms) : "                <<tms<<std::endl;

    // Utilities
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Utilities : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"tms.wait(2).count() : "                                                             <<tms.wait(2).count()<<std::endl;
    std::cout<<std::setw(width*2)<<"tms.wait(std::chrono::seconds(1)).count() : "                                       <<tms.wait(std::chrono::seconds(1)).count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.benchmark(10000000, [](){return 42.;}).count() : "                               <<ts.benchmark(1000000, [](){return 42.;}).count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.benchmark(10000000, [](double x){return std::sin(x);}, 42.).count() : "          <<ts.benchmark(1000000, [](double x){return std::sin(x);}, 42.).count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.benchmark(10000000, [](double x){return std::sin(x);}, 42).count() : "           <<ts.benchmark(1000000, [](double x){return std::sin(x);}, 42).count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.benchmark(10000000, [](double x){return std::lgamma(x);}, 42.).count() : "       <<ts.benchmark(1000000, [](double x){return std::lgamma(x);}, 42.).count()<<std::endl;
    std::cout<<std::setw(width*2)<<"ts.benchmark(10000000, std::plus<double>(), 4., 8.).count() : "                     <<ts.benchmark(1000000, std::plus<double>(), 4., 8.).count()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Timer::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // TIMER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
