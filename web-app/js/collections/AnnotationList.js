define(['backbone', 'models/Annotation'], 
    
        function(Backbone, Annotation) {
  'use strict';

  return Backbone.Collection.extend({
    model   : Annotation,
    url     : '/api/v1/annotations',

    getRange : function() {
      var range = []
      this.each(function(annotation) { range.push( _.range(annotation.get('start'), annotation.get('stop')+1)  ); })
      return _.uniq( _.flatten(range) );
    }

  });
});
