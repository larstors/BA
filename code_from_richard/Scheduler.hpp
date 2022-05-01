#ifndef SCHEDULER_HPP
#define SCHEDULER_HPP

#include <functional>
#include <utility>
#include <queue>
#include <vector>

// The Scheduler manages the waiting-time algorithm. Events are modelled as
// C++ function objects; these are called at the scheduled (simulation) time.
// The advance method advances to the next scheduled event and executes it.
class Scheduler {

public:

  using Event = std::function<void()>;
  using EventId = unsigned;

private:
  double _time = 0.0;

  using TimedEvent = std::pair<double, Event>;

  static constexpr auto cmp = [](const TimedEvent& a, const TimedEvent& b) {
    return a.first > b.first;
  };

  std::priority_queue<TimedEvent, std::vector<TimedEvent>, decltype(cmp)> queue;

public:

  Scheduler() : queue(cmp) { }

  // Schedule the event e to fire after the specified delay
  void schedule(double delay, Event e) {
    queue.push(std::make_pair(_time+delay, e));
  }

  // Advance the clock to the next event and execute it.
  // Returns true if an event is actually performed, otherwise false
  // (in this implementation, returns false only when the event queue is exhausted)
  bool advance() {
    if(queue.empty()) return false;
    auto te = queue.top();
    queue.pop();
    _time = te.first;
    te.second();
    return true;
  }

  // Advance the clock to the earlier of tmax and the next event time (if any).
  // If the next event falls before or at tmax, it will be executed.
  // Returns true if an event is actually performed, otherwise false.
  bool advance(double tmax) {
    if(tmax < _time) return false;
    if(queue.empty()) {
      _time = tmax;
      return false;
    }
    auto te = queue.top();
    if(te.first > tmax) {
      _time = tmax;
      return false;
    }
    queue.pop();
    _time = te.first;
    te.second();
    return true;    
  }

  auto time() const {
    return _time;
  }

  // Number of events pending
  auto pending() const {
    return queue.size();
  }

};

#endif
