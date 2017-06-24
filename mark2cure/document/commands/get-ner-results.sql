SELECT  `entity_recognition_entityrecognitionannotation`.`id`,
        `entity_recognition_entityrecognitionannotation`.`type_idx`,
        `entity_recognition_entityrecognitionannotation`.`text`,
        `entity_recognition_entityrecognitionannotation`.`start`,
        `document_annotation`.`created`,
        `document_document`.`id` as `document_pk`,
        `document_document`.`document_id` as `pmid`,
        `document_view`.`section_id`,
        `document_view`.`user_id`

FROM `entity_recognition_entityrecognitionannotation`

INNER JOIN `document_annotation`
    ON `document_annotation`.`object_id` = `entity_recognition_entityrecognitionannotation`.`id`
      AND `document_annotation`.`content_type_id` = {content_type_pk}

INNER JOIN `document_view`
    ON `document_annotation`.`view_id` = `document_view`.`id`

INNER JOIN `document_section`
    ON `document_view`.`section_id` = `document_section`.`id`

INNER JOIN `document_document`
    ON `document_document`.`id` = `document_section`.`document_id`

{filter_doc_level}

{filter_user_level}
