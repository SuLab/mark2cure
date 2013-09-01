define(['backbone', 'models/Annotation'], 

        function(Backbone, Annotation) {
  'use strict';

  return Backbone.Collection.extend({
    model   : Annotation,
    url     : '/api/v1/annotations',

    getRange : function() {
      //-- Returns back the indexes of the text which are part of a annotation
      var range = []
      this.each(function(annotation) { range.push( _.range(annotation.get('start'), annotation.get('stop')+1)  ); })
      return _.uniq( _.flatten(range) );
    },

    findContaining : function(index) {
      return this.filter(function(annotation) { 
        return index >= annotation.get('start') && index <= annotation.get('stop');
      });
    }

  });
});
