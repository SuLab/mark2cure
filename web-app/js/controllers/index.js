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

    showDocRelationship : function(doc_id) {
      vent.trigger('navigate:document:relationship', {doc_id: doc_id});
    },

    showDocument : function(doc_id, assignment_id, hit_id, worker_id, turk_sub) {
      // console.log('showDocument :: ', doc_id, assignment_id, hit_id, worker_id, turk_sub);
      User.set('assignment_id', null);
      window.aws = null;

      switch(assignment_id) {
      case undefined:
        //-- Normal user asking for specific document
        vent.trigger('navigate:document', {doc_id: doc_id});
        break;
      case 'ASSIGNMENT_ID_NOT_AVAILABLE':
        //-- Preview mode
        window.aws = {};
        window.aws.assignment_id = assignment_id;
        User.set('assignment_id', assignment_id);
        vent.trigger('navigate:document', {doc_id: doc_id});
        break;
      default:
        window.aws = {};
        window.aws.turk_sub = turk_sub;
        window.aws.worker_id = worker_id;
        window.aws.hit_id = hit_id;
        window.aws.assignment_id = assignment_id;
        window.aws.document_id = doc_id;
        //-- If via AMT, get that user started if not auth'd already
        User.set('assignment_id', assignment_id);

        if( User.authenticated() && User.get('mturk') ) {
          vent.trigger('navigate:document', {doc_id: doc_id});
        } else {
          User.set('username', worker_id);
          User.set('mturk', true);
          User.save(null, {success: function() {
            //-- After our user is saved, go ahead to get the document
            vent.trigger('navigate:document', {doc_id: doc_id});
          }});
        }
        break;
      }

    },

  };
});
