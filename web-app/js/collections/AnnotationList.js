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
    },

    exactMatch : function(word) {
      return this.filter(function(annotation) {
        return  word.get('start') == annotation.get('start') &&
                word.get('stop') == annotation.get('stop');
      });
    },

    add : function(ann) {
      //-- Prevent duplicate annotations from being submitted
      var isDupe = this.any(function(_ann) {
          return _ann.get('text') === ann.get('text') && _ann.get('start') === ann.get('start');
      });
      if (isDupe) { return false; }
      Backbone.Collection.prototype.add.call(this, ann);
    }

  });
});
