define(['vent', 'models/User'],
    function (vent, User) {
  'use strict';

  return {
    showSettings : function() { vent.trigger('navigate:settings', {}); },
    showInstructions : function() { vent.trigger('navigate:instructions', {}); },
    showSurvey : function() { vent.trigger('navigate:survey', {}); },
    sendMessage : function() { vent.trigger('navigate:message', {}); },

    showLibrary : function(quest) {
      vent.trigger('navigate:library', {'quest': quest});
    },

    showDocument : function(doc_id, assignment_id, hit_id, worker_id, turk_sub ) {
      console.log('showDocument', doc_id, assignment_id, hit_id, worker_id, turk_sub );

      if(assignment_id==undefined) {
        //-- Normal user asking for specific document
        vent.trigger('navigate:document', {doc_id: doc_id});
      } else {
        //-- If via AMT, get that user started
        User.save();
      }

      if(assignment_id == 'ASSIGNMENT_ID_NOT_AVAILABLE') {
        //-- Preview mode
        vent.trigger('navigate:document', {doc_id: doc_id});

      } else {
        vent.trigger('navigate:document', {doc_id: doc_id});

      }
    },

  };
});
