define(['marionette', 'bootstrap'], 
        function(marionette, bootstrap) {
  'use strict'

  return marionette.Region.extend({
    el: '#modal',

    constructor: function(){
      _.bindAll(this);
      marionette.Region.prototype.constructor.apply(this, arguments);
      //--  listens to the region’s “view”show” event, which is fired when a view’s contents are populated in to the `el` of the Region
      this.on('view:show', this.showModal, this);
    },

    getEl: function(selector){
      var $el = $(selector);
      $el.on('hidden', this.close);
      return $el;
    },

    showModal: function(view) {
      if (view.options.swap != true) {
        view.on('close', this.hideModal, this);
      };

      var type = true;
      if( view.options.static ) {
        type = 'static';
      }

      this.$el.modal({
        show : true,
        backdrop : type
      });
    },

    hideModal: function(){
      this.$el.modal('hide');
    },

  });
});
