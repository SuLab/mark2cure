define(['backbone', 'vent',
        //-- Data
        'models/Word', 'collections/WordList',
        'models/Annotation', 'collections/AnnotationList',
        'models/User',
        //-- Utilities
        'underscore.string', 'tastypie'],
    function(Backbone, vent,
        //-- Data
        Word, WordList,
        Annotation, AnnotationList,
        User,
        //-- Utilities
        _s){
  'use strict';

  //
  //-- This model repersents the concept of a Document
  //-- -- It stores the title, original text, id and any search metadata
  //

  return Backbone.RelationalModel.extend({
    url     : function() {
      return '/api/v1/documents/' + this.id;
    },
    defaults: {
      title         : '',
      text          : '',
      document_id   : 0,

      //-- for now, 'complete' can live entirely clientside
      complete      : false,

      matches       : []
    },

    relations: [{
      type: 'HasMany',
      key: 'annotations',

      relatedModel: Annotation,
      collectionType: AnnotationList,

      reverseRelation : {
        key : 'parentDocument',
        includeInJSON: true,
      }
    }, {
      type: 'HasMany',
      key: 'words',

      relatedModel: Word,
      collectionType: WordList,

      reverseRelation : {
        key : 'parentDocument',
        includeInJSON: false,
      }
    }],

    initialize : function() {
      this.bind('change:complete', this.evaluateResults);
      this.bind('change:complete', this.checkProgress);
    },

    parseText : function() {
      var self = this,
          step = 0,
          length = 0,
          words = _.map(_s.words( self.get('text') ), function(word) {
            length = word.length;
            step = step + length + 1;
            return {
              'text'      : word,
              'length'    : length,
              'start'     : step - length - 1,
              'stop'      : step - 2
            }
          });

        _.each(self.get('words'), function(word) {
          word.destroy();
        });
        self.get('words').add(words);
    },

    checkProgress : function() {
      if( this.collection && this.collection.completed().length == 5 ) {
        vent.trigger('modal:show_complete', {});
      }

      // this.collection.completed().length()/this.collection.length > .5
      // this.collection.fetchNext()

    },

    evaluateResults : function() {
      var self = this;

      //-- If the model was just completed
      if(this.get('complete') && this.hasChanged('complete') ) {
        $.getJSON('/api/v1/gm/'+this.id, function(data) {
          var gm = self.mapAnnotationsForComparision(data.objects);
          var user = self.mapAnnotationsForComparision(self.get('annotations').toJSON());

          var matches = _.intersectionObjects(gm, user);
          if(gm.length && user.length && matches.length) {
            self.set('matches', matches);
          }
        });

      };
    },

    mapAnnotationsForComparision : function(arr) {
      return _.map(arr, function(model) {
        return {'text'  : model.text,
                'start' : model.start}
          })
    },

    aws : function() {
      var self = this;
      console.log('aws :: ', User);
      if(User.authenticated() && User.get('mturk') && User.get('assignment_id')) {
        //-- Tell amazon they completed the hit
        $.ajax({
          type: "POST",
          url: "https://www.mturk.com/mturk/externalSubmit",
          data: { assignmentId : User.get('assignment_id'),
                  document_id : self.id,
                  annotations: self.get('annotations').pluck('text').join(', ')}
        });
      }

    }

  });
});
