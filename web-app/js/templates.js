define(function(require) {
  'use strict';
  return {
    header : {
      index      : require('tpl!templates/header/index.tmpl'),
    },

    analytics : {
      network : require('tpl!templates/analytics/network.tmpl')
    },

    main : {
      library : {
        index   : require('tpl!templates/main/library/index.tmpl'),
        list   : require('tpl!templates/main/library/list.tmpl'),
        item    : require('tpl!templates/main/library/item.tmpl')
      },

      game : {
        index : require('tpl!templates/main/game/index.tmpl'),
        controls : require('tpl!templates/main/game/controls/index.tmpl'),

        entity_tag : {
          worditem  : require('tpl!templates/main/game/entity_tag/worditem.tmpl')
        },

        results : {
          index     : require('tpl!templates/main/game/results/index.tmpl'),
        }
      }

    },

    footer : {
      index      : require('tpl!templates/footer/index.tmpl'),
    },

    modals : {
      first_run   : require('tpl!templates/modals/first_run.tmpl'),
      complete    : require('tpl!templates/modals/complete.tmpl'),

      survey      : require('tpl!templates/modals/survey.tmpl'),
      message     : require('tpl!templates/modals/message.tmpl'),

      settings      : require('tpl!templates/modals/settings.tmpl'),
      instructions  : require('tpl!templates/modals/instructions.tmpl')
    },

    snippets : {
      paragraph_info  : require('tpl!templates/snippets/paragraph_info.tmpl'),
      network_info    : require('tpl!templates/snippets/network_info.tmpl')
    }

  };
});
