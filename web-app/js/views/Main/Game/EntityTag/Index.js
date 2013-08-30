define(['marionette', 'templates', 'vent',
        //-- Views
        'views/Main/Game/EntityTag/Words'],
    function (Marionette, templates, vent,
              //-- Views
              Words) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.main.game.entity_tag.index,

    regions : {
      text   : 'p.paragraph'
    },

    onRender : function() {
      this.text.show( new Words({collection: this.model.get('words')}) );
    },

    //
    //-- Events
    //
    // saveGame : function() {
    //   //-- Annotations where sync'd with the server in real time
    //   _.each(this.model.get('words').getSelected(), function(word) {
    //     self.model.get('annotations').create({
    //       kind    : 0,
    //       type    : 'disease',
    //       position  : word.get('position'),
    //       text      : word.get('text'),
    //       length    : word.get('length'),
    //       start     : word.get('start'),
    //       stop      : word.get('stop')
    //     })
    //   });
    // }

  });
});
