define(['backbone', 
        //-- Models & Collections
        'models/Document'],

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

    parse : function(resp, xhr) {
      return resp.objects;
    },

    completed : function() {
      // return this.reset( this.where({complete : true}) );
      return this.where({complete : true});
    }

  });
});
