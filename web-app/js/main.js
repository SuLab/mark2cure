require.config({
  paths : {
    //-- jQuery Core
    'jquery'            : 'lib/jquery-1.9.1',

    //-- System Core
    'underscore'        : 'lib/underscore',
    'underscore.string' : 'lib/underscore.string',

    'backbone'          : 'lib/backbone',
    'relational'        : 'lib/backbone-relational',
    'marionette'        : 'lib/backbone.marionette',
    'backbone.flask'    : 'lib/backbone.flask',

    //-- UI Core
    'bootstrap'         : 'lib/bootstrap',
    'wizard'            : 'lib/wizard',
    'tourist'           : 'lib/tourist',

    //-- Utils
    'tpl'               : 'lib/tpl',
    'cookie'            : 'lib/jquery.cookie',
    'moment'            : 'lib/moment',
    'hash'              : 'lib/md5',
    'd3'                : 'lib/d3.v3',
    'b64'               : 'lib/webtoolkit.base64'
  },
  shim : {
    'd3' : {
      exports : 'd3'
    },

    'underscore' : {
      exports : '_'
    },

    'underscore.string' : {
        deps: ['underscore'],
        exports: '_s'
    },

    'relational' : {
      exports : 'Relational',
      deps: ['backbone', 'underscore']
    },

    'backbone' : {
      exports : 'Backbone',
      deps : ['jquery', 'underscore']
    },

    'marionette' : {
      exports : 'Backbone.Marionette',
      deps : ['backbone', 'relational']
    },

    'bootstrap' : {
      exports : 'bootstrap',
      deps    : ['jquery']
    },

  },
  deps : ['jquery','underscore', 'underscore.string']
});

require(
    ['app', 'backbone', 'routers/index', 'controllers/index'],
    function(app, Backbone, Router, Controller) {
  'use strict';
  app.start();
  new Router({
    controller : Controller
  });
  Backbone.history.start();
});
