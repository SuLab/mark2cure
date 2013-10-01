define(['backbone', 'relational'],
    function(Backbone){
  'use strict';

  return Backbone.RelationalModel.extend({
    defaults: {
      text      : '',

      length  : null,
      start   : null,
      stop    : null,

      kind    : 0, //-- entities
                    //--   [TX, "Type", [[START, STOP]] ]
                    //-- attributes
                    //--   [AX, "Type", TX]
                    //-- relations
                    //--   [ RX, "Type", [["Text", TX], ["Text", TX]] ]
                    //-- triggers
                    //--
                    //-- events
      type    : '',
    },

    sync    : function () { return false; },

  });
});
