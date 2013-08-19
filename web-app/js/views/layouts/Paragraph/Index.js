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
      this.preRender();

      _.bindAll(this, 'keyPressed', 'keyReleased');
      $(document).bind('keydown', this.keyPressed);
      $(document).bind('keyup', this.keyReleased);
      // $(document).bind('mouseup', this.mouseReleased);
    },

    onRender : function() {
      this.text.show( new Words(this.options) );
    },

    onBeforeClose : function() {
     // $(document).unbind('keypress', 'on_keypress');
    },

    preRender : function() {
      var self = this,
          //------
          // HEATMAP CUTOFF %
          // -----
          community_threshold = .65;
      this.options.max_pop = _.max( this.model.get('popularity') );
      this.options.color_scale = d3.scale.linear()
                                  .domain([0, this.options.max_pop])
                                  .range(["white", "yellow"]);

      if( this.model.get('complete') ) {
        var defined_correct = _.map(this.model.get('popularity'), function(i) { return i/self.options.max_pop >= community_threshold ? true : false; }),
            selected_index  = this.model.get('annotations').pluck('position'),
            user_correct    = _.map(this.model.get('popularity'), function(x, index) { return _.contains(selected_index, index); });
        this.options.score = this.compare(defined_correct, user_correct);
        this.options.ann_range = this.model.get('annotations').getRange();
      }
    },

    //
    //-- Events
    //
    keyPressed : function(evt) {
      if( evt.keyCode === 17 ) { this.model.get('cache').set('keypress', true); };
    },

    keyReleased : function(evt) {
      if( evt.keyCode === 17 ) { this.model.get('cache').set('keypress', false); };
    },

    mouseReleased : function(evt) {
      if(!document.getSelection) { return false; }
      var sel = document.getSelection(),
          $sel = $(sel.focusNode.parentElement),
          sel_text = sel.toString();

      if( $sel.closest('p').attr('class') === "paragraph editing" &&
          sel_text.length ) {
      }
    },

    submitAnnotations : function(evt) {
      var self = this;
      evt.preventDefault();
      this.saveAnnotations();
      this.render();

      if( this.collection.completed().length == 5 ) {
        vent.trigger('navigate:show_complete');
      };

      if( this.collection.completed().length == 1 ) {
        // vent.trigger('navigate:explain stuff here');
        this.ui.navigate.popover({  title   : 'Congrats! You annotated your first document!',
                                    content : '<p>Do your annotations agree with other users? Yellow highlights show community consensus</p><button class="btn btn-primary btn-small explain-network" type="button">Next</button>',
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
      this.ui.navigate.popover('hide');
      $('button.close').popover('hide');

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

      // if( $('#network').length ) {
      //     vent.trigger('navigate:analytics', {toggle: true});
      // }
    },

    //
    //-- Utils
    //
    compare : function(user, truth) {
      var a = 0, b = 0, c = 0, d = 0, fpr, fnr;
      if( truth.length !== user.length ) return;
      for (var i = 0; i < truth.length; i++) {
        switch(user[i] + truth[i]) {
          case 0: ++d; break;
          case 2: ++a; break;
          case 1: user[i] ? ++c : ++b; break;
        }
      }
      if( a+b+c+d !== truth.length ) return;
      fpr = c/(a+c); fnr = b/(b+d);
      return {'true_pos': a, 'false_neg': b, 'false_pos': c, 'true_neg': d, 'false_pos_rate': fpr, 'false_neg_rate': fnr}
    }

  });
});
