// $Id: DiskWriterResources.h,v 1.5 2009/10/13 15:08:33 mommsen Exp $
/// @file: DiskWriterResources.h 


#ifndef EventFilter_StorageManager_DiskWriterResources_h
#define EventFilter_StorageManager_DiskWriterResources_h

#include "EventFilter/StorageManager/interface/ErrorStreamConfigurationInfo.h"
#include "EventFilter/StorageManager/interface/EventStreamConfigurationInfo.h"

#include "boost/thread/condition.hpp"
#include "boost/thread/mutex.hpp"

#include <deque>
#include <stdint.h>


namespace stor
{

  /**
   * Container class for resources that are needed by the DiskWriter
   * and need to be accessed from multiple threads.
   *
   * $Author: mommsen $
   * $Revision: 1.5 $
   * $Date: 2009/10/13 15:08:33 $
   */

  class DiskWriterResources
  {
  public:

    /**
     * Constructs a DiskWriterResources instance.
     */
    DiskWriterResources();

    /**
     * Destructor.
     */
    virtual ~DiskWriterResources() {}

    /**
     * Requests that the DiskWriter streams be configured with the
     * specified configurations.  Also allows a new dequeue timeout
     * value to be specified. Existing stream configurations will be
     * discarded. 
     */
    void requestStreamConfiguration
    (
      EvtStrConfigListPtr,
      ErrStrConfigListPtr,
      double timeoutValue
    );

    /**
     * Requests that the DiskWriter streams be destroyed.
     * Any pending stream configuraton request will be cancelled.
     */
    void requestStreamDestruction();

    /**
     * Checks if a request has been made to change the stream configuration
     * in the DiskWriter streams *and* clears any pending request.
     * Existing streams are purged.
     * If doConfig is true, the supplied new configurations and a new dequeue
     * timeout value should be used to configure a new set of DiskWriter streams.
     */
    bool streamChangeRequested
    (
      bool& doConfig,
      EvtStrConfigListPtr&,
      ErrStrConfigListPtr&,
      double& timeoutValue
    );

    /**
     * Waits until a requested stream configuration change has been completed.
     */
    virtual void waitForStreamChange();

    /**
     * Returns true when a stream configuration change is requested or in progress.
     */
    virtual bool streamChangeOngoing();

    /**
     * Indicates that the stream configuration change is done.
     */
    void streamChangeDone();

    /**
     * Sets the DiskWriter "busy" status to the specified state.
     */
    void setBusy(bool isBusyFlag);

    /**
     * Tests if the DiskWriter is currently busy..
     */
    bool isBusy();

  private:

    bool _configurationIsNeeded;
    bool _streamChangeIsNeeded;
    bool _fileClosingTestIsNeeded;
    bool _diskWriterIsBusy;

    EvtStrConfigListPtr _requestedEventStreamConfig;
    ErrStrConfigListPtr _requestedErrorStreamConfig;
    double _requestedTimeout;

    bool _streamChangeInProgress;
    boost::condition _streamChangeCondition;
    
    mutable boost::mutex _streamChangeMutex;
  };

}

#endif

/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
