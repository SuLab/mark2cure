define(['marionette', 'templates', 'vent',
        //-- Models & Collections
        'views/layouts/Paragraph/Words',
        //-- ETC
        'd3'], 

        function (Marionette, templates, vent,
                  Words) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.paragraph.index,
    templateHelpers : function() { return this.options; },

    className : 'text-view',

    regions : {
      text      : 'p.paragraph',
    },

    ui : {
      'complete_alert'  : '.alert.alert-success',
      'navlist'         : 'ul.pager',
      'navigate'        : 'button.navigate'
    },

    events : {
      'click button.done'   : 'submitAnnotations',
      'click button.navigate'    : 'nextDocument',
    },

    initialize : function(options) {
      //-- This is the Layout view to manage the Paragraph Views/Logic
      //-- this.model == The Document Resource
      //-- this.collection == The Collection of all Documents
      //-- options.user == The Currently Logged in User
      //--
      //-- options.max_pop == Highest ranked word
      //-- options.score   == Result score of user's selection
      options.max_pop = 0; options.score = {}; options.ann_range = [];
    },

    onRender : function() {
      this.text.show( new Words(this.options) );
    },

    //
    //-- Events
    //
    submitAnnotations : function(evt) {
      var self = this;
      evt.preventDefault();
      this.saveAnnotations();
      this.render();

      if( this.collection.completed().length == 5 ) {
        vent.trigger('navigate:show_complete');
      };

      if( this.collection.completed().length == 1 ) {
        this.ui.navigate.popover({  title   : 'Congrats! You annotated your first document!',
                                    content : templates.snippets.paragraph_info({}),
                                    html    : true,
                                    trigger : 'manual',
                                    placement : 'left',
                                    container : 'body' });
        this.ui.navigate.popover('show');

        $('button.explain-network').click(function(e) {
          self.ui.navigate.popover('hide');
          vent.trigger('navigate:analytics', {toggle: true, explain: true});
        });

        this.options.user.save({'first_run' : false});
      }

    },

    nextDocument : function(evt) {
      //-- Incase they didn't close it before or are skipping
      $('.popover').hide()

      evt.preventDefault();
      var index = this.collection.indexOf(this.model) + 1,
          index = (index < 0) ? this.collection.length-1 : index,
          index = (index == this.collection.length) ? 0 : index;

      this.options.model = this.model = this.collection.at(index);
      Backbone.history.navigate( '#/'+ this.model.id );
      this.preRender();
      this.render();
    },

    saveAnnotations : function() {
      var self = this;

      //-- Annotations where sync'd with the server in real time
      _.each(this.model.get('words').getSelected(), function(word) {
        self.model.get('annotations').create({
          kind    : 0,
          type    : 'disease',

          position  : word.get('position'),
          text      : word.get('text'),
          length    : word.get('length'),
          start     : word.get('start'),
          stop      : word.get('stop')
        })
      });

      this.model.save({'complete': true});
      this.preRender();
      this.render();

      if( $('#network').length ) { vent.trigger('network:refresh', {}); }
    },

  });
});
