define(['backbone', 
        //-- Models & Collections
        'models/Document', 'tastypie'],

        function( Backbone, 
                  //-- Models & Collections
                  Document) {
  'use strict';

  return Backbone.Collection.extend({
    model   : Document,
    url     : '/api/v1/documents',

    // comparator: function(doc) {
      // return -doc.get('created');
    // },

    completed : function() {
      // return this.reset( this.where({complete : true}) );
      return this.where({complete : true});
    },

  });
});
